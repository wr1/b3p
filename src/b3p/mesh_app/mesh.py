# Mesh building functions for mesh_app, self-contained but synced with build mesh interpolation.

import os
import logging
import pickle
from copy import deepcopy as dc
import numpy as np
from pathlib import Path
import pyvista as pv
import vtk
from ..geometry.blade import blade
from ..geometry.blade_section import section as GeometrySection
from ..geometry.loft_utils import load, interp, optspace
from ..geometry.splining import intp_c

logger = logging.getLogger(__name__)

# Local classes updated to match build mesh interpolation

def spline_interp(x, y, newx):
    spl = vtk.vtkCardinalSpline()
    spl.SetLeftConstraint(2)
    spl.SetRightConstraint(2)
    for i in zip(x, y):
        spl.AddPoint(i[0], i[1])
    return [spl.Evaluate(i) for i in newx]


def distance(point1, point2):
    return np.sqrt(sum((i[1] - i[0]) ** 2 for i in zip(point1, point2)))


def equals(v1, v2):
    tol = 1e-6
    return (v1 - v2) ** 2 < tol


def mesh_line(pnt1, pnt2, n_cells, id):
    xyz = []
    web_height = vtk.vtkGeoMath().DistanceSquared(pnt1, pnt2) ** 0.5
    for i in zip(pnt1, pnt2):
        mm = min(0.3, 0.06 / web_height)
        rel = sorted([0, 1] + list(np.linspace(mm, 1.0 - mm, n_cells - 2)))
        ab = [j * (i[1] - i[0]) + i[0] for j in rel]
        xyz.append(np.array(ab))
    dst = [i[1:] - i[:-1] for i in xyz]
    sl = (dst[0] ** 2 + dst[1] ** 2 + dst[2] ** 2) ** 0.5
    pl = [0] + [sum(sl[:i]) for i in range(1, len(sl) + 1)]
    ppl = [-i + pl[-1] for i in pl]
    ml = [abs(i - 0.5 * web_height) for i in pl]
    wh = [web_height for _ in ml]
    rad = np.mean(xyz[2])
    r = [rad for _ in range(n_cells)]
    arrays = {
        "d_te": pl,
        "d_le": ppl,
        "d_le_r": [i / max(ppl) for i in ppl],
        f"d_{id}_r": [i / max(ml) for i in ml],
        f"d_{id}": ml,
        "web_height": wh,
        "radius": r,
        "is_web": [1.0 for _ in ppl],
    }
    return list(zip(*xyz)), arrays


class GeometrySection:
    def __init__(self, r, r_relative, points, min_te_thickness=0.002, open_te=False):
        self.r = r
        self.r_relative = r_relative
        if open_te:
            self.base_points = self._open_te(points, min_te_thickness)
        else:
            self.base_points = points
        self.te_thickness = min_te_thickness
        self.points = vtk.vtkPoints()
        for i in self.base_points:
            self.points.InsertNextPoint((i[0], i[1], r))
        self.poly = vtk.vtkPolyData()
        self.poly.SetPoints(self.points)

    def _open_te(self, points, te_thickness):
        if len(points[0]) == 2:
            points = [list(i) + [0] for i in points]
        else:
            points = [list(i) for i in points]
        if (
            vtk.vtkMath.Distance2BetweenPoints(points[0], points[-1]) ** 0.5
            < te_thickness
        ):
            points[-1][1] -= (
                te_thickness
                - vtk.vtkMath.Distance2BetweenPoints(points[0], points[-1]) ** 0.5
            )
            for i in reversed(points):
                d = vtk.vtkMath.Distance2BetweenPoints(points[-1], i) ** 0.5
                if d >= 0.1:
                    break
                if i != points[-1]:
                    points.remove(i)
        return points

    def set_twist(self, twist):
        self.twist = twist
        self.poly.rotate_z(twist, inplace=True)

    def translate(self, dx, dy):
        self.poly.translate([dx, dy, 0.0], inplace=True)

    def scale(self, sx, sy):
        self.poly.scale([sx, sy, 1.0], inplace=True)

    def _create_evaluations(
        self, n_points, webs, mult=True, make_plots=False, panel_mesh_scale=None
    ):
        if panel_mesh_scale is None:
            panel_mesh_scale = []
        if webs == []:
            return np.linspace(0.0, 1.0, n_points)
        avspl = []
        for i in webs:
            avspl.extend(i.average_splits())
        intervals = sorted([0, 1] + list(avspl))
        isize = [intervals[i + 1] - intervals[i] for i in range(len(intervals) - 1)]
        for ps in panel_mesh_scale:
            if ps[0] < len(isize):
                isize[ps[0]] *= ps[1]
        inum = [int(round(n_points * i / sum(isize))) for i in isize]
        if mult:
            isize[0] *= 2
            isize[-1] *= 2
        for _ in range(5):
            if sum(inum) < n_points:
                inum[inum.index(min(inum))] += 1
            elif sum(inum) > n_points:
                inum[inum.index(max(inum))] -= 1
        splts = []
        for i in webs:
            splts.extend(i.splits(self.r, self.r_relative))
        real_intervals = sorted([0, 1] + list(splts))
        pnts = []
        for i in range(len(inum)):
            interval = np.linspace(
                real_intervals[i], real_intervals[i + 1], inum[i] + 1 * (i != 0)
            )
            pnts.extend(interval)
        return sorted(set(pnts))

    def respline(self, n_points, webs=None, added_datums=None, panel_mesh_scale=None):
        if webs is None:
            webs = []
        if added_datums is None:
            added_datums = {}
        if panel_mesh_scale is None:
            panel_mesh_scale = []
        # Updated to match build mesh interpolation: use vtkParametricSpline with constraints 3
        spline = vtk.vtkParametricSpline()
        spline.SetPoints(self.poly.GetPoints())
        spline.SetLeftConstraint(3)
        spline.SetLeftValue(1.0)
        spline.SetRightConstraint(3)
        spline.SetRightValue(1.0)
        p, du = [0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]
        spline.DerivativesAvailableOn()
        out = []
        last = []
        distance_along_airfoil = 0.0
        dist = []
        dist_from_te = []
        rel_dist_from_te = []
        x_abs, y_abs = [], []
        te = [
            (i[0] + i[1]) / 2.0
            for i in zip(
                self.poly.GetPoint(0),
                self.poly.GetPoint(self.poly.GetNumberOfPoints() - 1),
            )
        ]
        evaluations = self._create_evaluations(
            n_points, webs, panel_mesh_scale=panel_mesh_scale
        )
        for i in evaluations:
            p[0] = i
            u = [0, 0, 0]
            spline.Evaluate(p, u, du)
            out.append(u)
            if last == []:
                last = u
            distance_along_airfoil += distance(u, last)
            dist_from_te.append(distance(u, te))
            rel_dist_from_te.append(i)
            dist.append(distance_along_airfoil)
            x_abs.append(u[0])
            y_abs.append(u[1])
            last = u
        x, y, z = zip(*out)
        dist_miny = dist[y.index(min(y))]
        dist_sparcap = [i - dist_miny for i in dist]
        max_dist_from_te = max(dist_from_te)
        le_point = dist_from_te.index(max_dist_from_te)
        cline = [i[1] - i[0] for i in zip(out[le_point], te)]
        nor = [0, 0, 0]
        vtk.vtkMath.Cross(cline, [0, 0, -1], nor)
        chord = []
        for i in out:
            v1 = cline
            v2 = nor
            v3 = [j[1] - j[0] for j in zip(te, i)]
            chord.append(
                max(
                    0,
                    -(v3[0] / v1[0] - (v2[0] * v3[1]) / (v2[1] * v1[0]))
                    / (1.0 - (v2[0] * v1[1]) / (v2[1] * v1[0])),
                )
            )
        lechord = [-i + 1.0 for i in chord]
        lechord_absolute = [i * (cline[0] ** 2 + cline[1] ** 2) ** 0.5 for i in lechord]
        techord_abs = [i * (cline[0] ** 2 + cline[1] ** 2) ** 0.5 for i in chord]
        le_datum = [abs(dist[le_point] - i) for i in dist]
        te_datum = [
            (dist[i] if i < le_point else max(dist) - dist[i]) for i in range(len(dist))
        ]
        web_datums = {}
        for i in webs:
            splits = sorted(i.splits(self.r, self.r_relative))
            datum, datum_r = [], []
            for j in zip(rel_dist_from_te, [1 if k < len(rel_dist_from_te) // 2 else -1 for k in range(len(rel_dist_from_te))]):
                datum.append(
                    (j[0] - splits[j[1]]) * distance_along_airfoil * (-1 if j[1] else 1)
                )
                datum_r.append((j[0] - splits[j[1]]) * (-1 if j[1] else 1))
            web_datums[f"d_{i.coordinate}"] = datum
            web_datums[f"d_{i.coordinate}_r"] = datum_r
        r = [self.r for _ in dist]
        mdist = max(dist)
        datums = {
            "d_te": te_datum,
            "radius": r,
            "rr": self.r_relative * np.ones_like(r),
            "d_rel_dist_from_te": rel_dist_from_te,
            "d_abs_dist_from_te": dist,
            "d_abs_dist_from_bte": [-i + mdist for i in dist],
            "chord_length": [max_dist_from_te for _ in dist_from_te],
            "d_miny": dist_sparcap,
            "d_le": le_datum,
            "zone_ss": [1 if i < len(dist) // 2 else -1 for i in range(len(dist))],
            "zone_ps": [0 if i < len(dist) // 2 else 1 for i in range(len(dist))],
            "d_te_r": [i / max(te_datum) if max(te_datum) > 0 else 0 for i in te_datum],
            "d_le_r": [i / max(le_datum) if max(le_datum) > 0 else 0 for i in le_datum],
            "d_chord": chord,
            "d_techord": chord,
            "d_techord_abs": techord_abs,
            "d_lechord": lechord,
            "d_lechord_abs": lechord_absolute,
            "d_x": x_abs,
            "d_y": y_abs,
            "d_sl": [dist[i] - dist[len(dist) // 2] for i in range(len(dist))],
            "d_sla": [abs(dist[i] - dist[len(dist) // 2]) for i in range(len(dist))],
            "is_web": [0 for _ in r],
        }
        for i in web_datums:
            datums[i] = np.array(web_datums[i]).astype(np.float32)
        for i in added_datums.items():
            offs = np.interp(self.r_relative, i[1][1], i[1][2])
            datums[i[0]] = np.array(np.array(datums[i[1][0]]) + offs).astype(np.float32)
        return out, datums


class Web:
    def __init__(
        self, points, web_root, web_tip, web_name, coordinate, flip_normal=False
    ):
        self.points = points
        self.web_root = web_root
        self.web_tip = web_tip
        spl1, spl2 = vtk.vtkSCurveSpline(), vtk.vtkSCurveSpline()
        for i in points:
            spl1.AddPoint(i[0], i[1])
            spl2.AddPoint(i[0], i[2])
        self.splines = (spl1, spl2)
        self.evaluations = {}
        self.name = web_name
        self.coordinate = coordinate
        self.flip_normal = flip_normal

    def average_splits(self):
        g = list(zip(*self.points))
        return np.mean(g[1]), np.mean(g[2])

    def splits(self, r, r_relative):
        out = (self.splines[0].Evaluate(r), self.splines[1].Evaluate(r))
        self.evaluations[int(round(r * 1e2) * 10)] = [out]
        return out

    def _find_top_and_bottom_points(self, mesh):
        rad = mesh.GetPointData().GetArray("radius")
        rel_dist = mesh.GetPointData().GetArray("d_rel_dist_from_te")
        for i in range(mesh.GetNumberOfPoints()):
            rm = rad.GetValue(i)
            if self.web_root <= rm <= self.web_tip:
                rd = rel_dist.GetValue(i)
                pnt = mesh.GetPoint(i)
                rmm = int(round(rm * 1e2) * 10)
                if equals(rd, self.evaluations[rmm][0][0]) or equals(
                    rd, self.evaluations[rmm][0][1]
                ):
                    self.evaluations[rmm].append(pnt)

    def _create_quad_connectivity(self, n_points, n_total, flip=False):
        nrows = int(n_total / n_points)
        colids = np.arange(n_points - 1)
        npp = (
            np.arange(1, nrows).repeat(n_points - 1).reshape(nrows - 1, n_points - 1)
            - 1
        ) * n_points + colids
        if flip:
            stck = [
                np.ones_like(npp) * 4,
                npp + 1,
                npp + n_points + 1,
                npp + n_points,
                npp,
            ]
        else:
            stck = [
                np.ones_like(npp) * 4,
                npp,
                npp + n_points,
                npp + n_points + 1,
                npp + 1,
            ]
        return np.stack(stck).T.flatten()

    def _create_points(self, n_cells):
        ev = self.evaluations
        vp = []
        added_arrays = {}
        for i in sorted(ev):
            if len(ev[i]) == 3:
                pnts, data = mesh_line(ev[i][1], ev[i][2], n_cells, self.coordinate)
                vp.extend(pnts)
                for j in data:
                    if j not in added_arrays:
                        added_arrays[j] = data[j]
                    else:
                        added_arrays[j].extend(data[j])
        return vp, added_arrays

    def write_mesh(self, vtpfile):
        if hasattr(self, "webmesh"):
            self.webmesh.save(vtpfile)
        logger.info(f"Wrote mesh to {vtpfile}")

    def mesh(self, mesh, n_cells):
        self._find_top_and_bottom_points(mesh)
        points, pdata = self._create_points(n_cells)
        if not points:
            self.webmesh = pv.PolyData()
            logger.info(f"Web {self.name}: no points, creating empty mesh")
            return
        cells = self._create_quad_connectivity(n_cells, len(points), self.flip_normal)
        self.webmesh = pv.PolyData(points, faces=cells)
        for i in pdata:
            self.webmesh.point_data[i] = np.array(pdata[i]).astype(np.float32)


class BladeShape:
    def __init__(
        self,
        sections=None,
        section_resolution=200,
        web_resolution=20,
        added_datums=None,
        prefix="",
    ):
        if sections is None:
            sections = []
        if added_datums is None:
            added_datums = {}
        self.set_sections(sections)
        self.set_section_resolution(section_resolution)
        self.webs = []
        self.web_resolution = web_resolution
        self.added_datums = added_datums
        self.prefix = prefix

    def set_web(self, web):
        self.webs.append(web)

    def set_sections(self, sections):
        self.sections = sections

    def set_section_resolution(self, n_points):
        self.n_points = n_points

    def build_interpolated_sections(self, radii, interpolation_type=1):
        r = [i.r for i in self.sections]
        pnts = [i.respline(self.n_points)[0] for i in self.sections]
        nxyz = []
        for i in zip(*pnts):
            if interpolation_type == 1:
                nx = spline_interp(r, list(zip(*i))[0], radii)
                ny = spline_interp(r, list(zip(*i))[1], radii)
                nz = spline_interp(r, list(zip(*i))[2], radii)
            elif interpolation_type == 2:
                nx = np.interp(radii, r, list(zip(*i))[0])
                ny = np.interp(radii, r, list(zip(*i))[1])
                nz = np.interp(radii, r, list(zip(*i))[2])
            elif interpolation_type == 3:
                nx = spline_interp_k(r, list(zip(*i))[0], radii)
                ny = spline_interp_k(r, list(zip(*i))[1], radii)
                nz = spline_interp_k(r, list(zip(*i))[2], radii)
            nxyz.append(list(zip(nx, ny, nz)))
        self.interp_sections = [
            GeometrySection(i[0], i[0] / max(radii), i[1])
            for i in zip(radii, zip(*nxyz))
        ]

    def mesh(self, n_points=100, close=True, panel_mesh_scale=None):
        if panel_mesh_scale is None:
            panel_mesh_scale = []
        self.poly = vtk.vtkPolyData()
        points = vtk.vtkPoints()
        added_arrays = {}
        for i in self.interp_sections:
            pnts, data = i.respline(
                n_points, self.webs, self.added_datums, panel_mesh_scale
            )
            for j in pnts:
                points.InsertNextPoint(j[0], j[1], j[2])
            for j in data:
                if j not in added_arrays:
                    added_arrays[j] = vtk.vtkFloatArray()
                    added_arrays[j].SetName(j)
                for k in data[j]:
                    added_arrays[j].InsertNextValue(k)
        self.poly.SetPoints(points)
        np_total = self.poly.GetNumberOfPoints()
        quads = vtk.vtkCellArray()
        for i in range(1, int(np_total / n_points)):
            np_start = (i - 1) * n_points
            nc_start = i * n_points
            for j in range(n_points - (0 if close else 1)):
                quads.InsertNextCell(4)
                quads.InsertCellPoint(np_start + j)
                quads.InsertCellPoint(nc_start + j)
                quads.InsertCellPoint(nc_start + (j + 1) % n_points)
                quads.InsertCellPoint(np_start + (j + 1) % n_points)
        self.poly.SetPolys(quads)
        for value in added_arrays.values():
            self.poly.GetPointData().AddArray(value)
        for i in self.webs:
            i.mesh(self.poly, self.web_resolution)

    def write_mesh(self, filename):
        try:
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName(filename)
            writer.SetInputData(self.poly)
            writer.Write()
        except Exception:
            logger.info("no valid mesh available")
        for i in self.webs:
            if hasattr(i, "webmesh"):
                i.write_mesh(f"{i.name}.vtp")


def build_blade_mesh(config, workdir):
    """Build the blade mesh using self-contained logic in mesh_app."""
    prefix = config["general"]["prefix"]
    radii = np.linspace(0, 100, 100)
    web_inputs = config["mesh"]["webs"]
    base_vtp = workdir / f"{prefix}_base.vtp"
    web_intersections = build_webs(
        str(base_vtp), web_inputs, prefix=prefix, workdir=workdir
    )
    pckfile = workdir / f"{prefix}.pck"
    outfile = workdir / "blade_mesh.vtp"
    build_mesh(
        str(pckfile),
        radii,
        web_inputs,
        web_intersections,
        prefix,
        outfile=str(outfile),
    )
    web_files = [f"{prefix}_{i}.vtp" for i in web_inputs]
    logger.info(f"Web planes exported to {workdir}: {', '.join(web_files)}")


def build_mesh(
    pckfile,
    radii,
    web_inputs,
    web_intersections,
    prefix,
    n_web_points=10,
    n_ch_points=120,
    outfile="out.vtp",
    added_datums=None,
    panel_mesh_scale=None,
):
    if added_datums is None:
        added_datums = {}
    if panel_mesh_scale is None:
        panel_mesh_scale = []
    sections = pickle.load(open(pckfile, "rb"))
    weblist = [
        Web(
            points=web_intersections[i],
            web_root=web_inputs[i]["z_start"],
            web_tip=web_inputs[i]["z_end"],
            web_name=f"{prefix}_{i}",
            coordinate=i,
            flip_normal=(web_inputs[i]["origin"][1] > 0),
        )
        for i in web_inputs
    ]
    nsec = []
    z = [i[0][2] for i in sections]
    for i in sections:
        r = i[0][2] - min(z)
        r_rel = (r - min(z)) / (max(z) - min(z))
        sec = GeometrySection(r, r_rel, i, open_te=False)
        nsec.append(sec)
    blade = BladeShape(
        nsec,
        section_resolution=200,
        web_resolution=n_web_points,
        added_datums=added_datums,
        prefix=prefix,
    )
    for i in weblist:
        blade.set_web(i)
    blade.build_interpolated_sections(radii=radii, interpolation_type=2)
    blade.mesh(n_ch_points, panel_mesh_scale=panel_mesh_scale)
    blade.write_mesh(outfile)
    logger.info(f"Wrote blade mesh to {outfile}")
    return blade


def build_webs(mesh_path, webs, prefix="__dum", workdir=Path(".")):
    """Build web intersections using pyvista and numpy with matrix operations, using radius array for reference."""
    mesh = pv.read(mesh_path)
    web_meshes = {}
    for i in webs:
        normal = (0, 1, 0)
        name = str(prefix) + "_" + i
        if "orientation" in webs[i]:
            normal = webs[i]["orientation"]
        loc = np.array(webs[i]["origin"])
        # Slice the mesh with the plane
        slice_mesh = mesh.slice(normal=normal, origin=loc)
        if slice_mesh.n_points == 0:
            logger.warning(f"Web {name}: no intersection points")
            continue
        # Get point data
        c = slice_mesh.point_data.get("d_rel_dist_from_te")
        cc = slice_mesh.point_data.get("d_abs_dist_from_te")
        p = slice_mesh.point_data.get("d_rel_dist_from_te")
        radius_coords = slice_mesh.point_data.get("radius")
        if c is None or cc is None or radius_coords is None:
            logger.warning(f"Web {name}: missing point data arrays")
            continue
        # Use radius_coords instead of z_coords for interpolation reference
        minr = np.min(radius_coords)
        maxr = np.max(radius_coords)
        r = np.linspace(minr, maxr, 400)
        # Use numpy boolean indexing for vectorized operations
        leading_mask = c > 0.5
        trailing_mask = ~leading_mask
        # Leading edge
        if np.any(leading_mask):
            lw_r = radius_coords[leading_mask]
            lw_p = p[leading_mask]
            lw1 = np.interp(r, lw_r, lw_p)
        else:
            lw1 = np.zeros_like(r)
        # Trailing edge
        if np.any(trailing_mask):
            ww_r = radius_coords[trailing_mask]
            ww_p = p[trailing_mask]
            ww1 = np.interp(r, ww_r, ww_p)
        else:
            ww1 = np.zeros_like(r)
        # Ratio
        rt1 = np.interp(r, radius_coords, cc / c)
        # Build out_list using vectorized operations
        z_follow_blade = webs[i]["z_follow_blade"]
        mask_follow = r <= z_follow_blade
        lwl = np.where(mask_follow, lw1, 0)
        wwl = np.where(mask_follow, ww1, 0)
        out_list = np.column_stack((r, wwl, lwl, rt1)).tolist()
        # Append end point if needed
        if webs[i]["z_end"] > out_list[-1][0]:
            out_list.append([webs[i]["z_end"], out_list[-1][1], out_list[-1][2], out_list[-1][3]])
        # Collect points for JSON and VTP
        lwp = slice_mesh.points[leading_mask].tolist()
        wwp = slice_mesh.points[trailing_mask].tolist()
        data = {
            "name": name,
            "data": out_list,
            "points": {"lwp": lwp, "wwp": wwp},
        }
        with open(workdir / f"{name}.json", "w") as f:
            import json
            json.dump(data, f, indent=4)
        # Generate VTP file
        web_points = np.array(lwp + wwp)
        if web_points.size > 0:
            web_poly = pv.PolyData(web_points)
            web_poly.save(workdir / f"{name}.vtp")
            logger.info(f"Wrote web VTP to {workdir / f'{name}.vtp'}")
        web_meshes[i] = out_list
    return web_meshes
