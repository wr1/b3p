# Mesh building functions for mesh_app.

import os
import logging
import pickle
from copy import deepcopy as dc
import numpy as np
from pathlib import Path
import pyvista as pv
import vtk
from ..geometry.blade import blade
from ..geometry.blade_section import section
from ..geometry.loft_utils import load, interp, optspace
from ..geometry.splining import intp_c

logger = logging.getLogger(__name__)


def write_web(
    loc,
    normal,
    mesh,
    name,
    rootcut=0.0,
    tipcut=100.0,
    tip=0.0,
    zone="d_rel_dist_from_te",
    workdir=Path("."),
):
    """Slices a plane through a mesh, and reports back a relative
    coordinate of top and bottom lines. This is used to represent geometrically straight entities
    in 3D as coordinates in local systems defined per section.

    Args:
        loc (np.ndarray): location of the plane
        normal (tuple): normal of the plane
        mesh (str): mesh to slice
        name (str): name of the web
        rootcut (float, optional): cut the root of the web. Defaults to 0.0.
        tipcut (float, optional): cut the tip of the web. Defaults to 100.0.
        tip (float, optional): tip of the blade. Defaults to 0.0.
        zone (str, optional): zone of the blade. Defaults to "d_rel_dist_from_te".
        workdir (Path, optional): work directory. Defaults to Path(".").

    Returns:
        list: A list of tuples, each containing (rr, wwl, lwl, rt1[i])
    """

    rd = vtk.vtkXMLPolyDataReader()
    rd.SetFileName(mesh)
    rd.Update()
    poly = rd.GetOutput()

    plane = vtk.vtkPlane()
    plane.SetOrigin(loc)
    plane.SetNormal(normal)

    clip = vtk.vtkCutter()
    clip.SetCutFunction(plane)
    clip.SetInputData(poly)

    clip.Update()
    out = clip.GetOutput()
    points = []

    c = out.GetPointData().GetArray("d_rel_dist_from_te")
    cc = out.GetPointData().GetArray("d_abs_dist_from_te")
    p = out.GetPointData().GetArray(zone)

    num_points = out.GetNumberOfPoints()
    logger.info(
        f"Web {name}: origin {loc}, normal {normal}, points found: {num_points}"
    )

    if c is None or p is None or cc is None:
        logger.warning(f"Web {name}: missing point data arrays")
        return []

    lw, ww = [], []
    lwp, wwp = [], []
    rta = []
    minr, maxr = 1000, -1000
    for i in range(out.GetNumberOfPoints()):
        if c.GetValue(i) > 0.5:
            lw.append((out.GetPoint(i)[2], p.GetValue(i)))
            lwp.append(out.GetPoint(i))
            rta.append((out.GetPoint(i)[2], cc.GetValue(i) / c.GetValue(i)))
        else:
            ww.append((out.GetPoint(i)[2], p.GetValue(i)))
            wwp.append(out.GetPoint(i))
        minr = min(minr, out.GetPoint(i)[2])
        maxr = max(maxr, out.GetPoint(i)[2])
        points.append(out.GetPoint(i))

    if not lw or not ww:
        logger.warning(f"Web {name}: no points in leading or trailing edge")
        return []

    lw = list(zip(*lw))
    ww = list(zip(*ww))
    rt = list(zip(*rta))
    # interpolate the results
    r = np.linspace(minr, maxr, 400)
    lw1 = np.interp(r, lw[0], lw[1])
    ww1 = np.interp(r, ww[0], ww[1])
    rt1 = np.interp(r, rt[0], rt[1])

    out = []
    for i, rr in enumerate(r):
        if rr <= tipcut:
            lwl, wwl = lw1[i], ww1[i]
        out.append((rr, wwl, lwl, rt1[i]))
        if i > 0 and r[i - 1] < rootcut and rr > rootcut:
            for j in range(i):
                out[j][1] = wwl
                out[j][2] = lwl

    out = sorted(out)
    if tip > out[-1][0]:
        out.append((tip, out[-1][1], out[-1][2], out[-1][3]))

        # Prepare data for JSON serialization
        data = {
            "name": name,
            "data": out,
            "points": {"lwp": [list(p) for p in lwp], "wwp": [list(p) for p in wwp]},
        }

        # Write data to a JSON file
        with open(workdir / f"{name}.json", "w") as f:
            import json

            json.dump(data, f, indent=4)

    # open("%s.txt" % name, "wb").write(str(out).encode("utf-8"))
    # open("%s_points.txt" % name, "wb").write(str([lwp, wwp]).encode("utf-8"))

    return out


def build_webs(mesh, webs, prefix="__dum", workdir=Path(".")):
    """Builds web meshes based on provided web definitions.
    Args:
        mesh (object): The base mesh object to which the webs will be attached.
        webs (dict): A dictionary defining the webs to be created.
            Each key in the dictionary represents the name of a web, and the
            corresponding value is a dictionary containing the web's properties,
            including:
                "origin" (list/array-like): The origin point of the web.
                "z_start" (float): The z-coordinate where the web starts.
                "z_follow_blade" (float): The z-coordinate where the web follows the blade.
                "z_end" (float): The z-coordinate where the web ends.
                "orientation" (list/array-like, optional): The normal vector
                    defining the orientation of the web. Defaults to (0, 1, 0).
        prefix (str, optional): A prefix to be added to the name of each web mesh.
            Defaults to "__dum".
        workdir (Path, optional): work directory. Defaults to Path(".").
    Returns:
        dict: A dictionary containing the generated web meshes. The keys of the
            dictionary are the names of the webs, and the values are the
            corresponding mesh objects.
    """

    web_meshes = {}
    for i in webs:
        normal = (0, 1, 0)
        name = str(prefix) + "_" + i

        if "orientation" in webs[i]:
            normal = webs[i]["orientation"]

        fea_web = write_web(
            np.array(webs[i]["origin"]),
            normal,
            mesh,
            name,
            rootcut=webs[i]["z_start"],
            tipcut=webs[i]["z_follow_blade"],
            tip=webs[i]["z_end"],
            workdir=workdir,
        )
        web_meshes[i] = fea_web

    return web_meshes


def build_blade_mesh(config, workdir):
    """Build the blade mesh including webs."""
    pln = config["planform"]
    radii = np.linspace(0, 100, 100)
    web_inputs = config["mesh"]["webs"]
    base_vtp = workdir / "blade_geometry.vtp"
    web_intersections = build_webs(
        str(base_vtp), web_inputs, prefix="blade", workdir=workdir
    )
    prefix = "blade"
    pckfile = workdir / "blade_geometry.pck"
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
    """Build the 3D mesh with webs linked to the shell."""
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


# Local classes for self-containment


def spline_interp(x, y, newx):
    spl = vtk.vtkCardinalSpline()
    spl.SetLeftConstraint(2)
    spl.SetRightConstraint(2)
    for i in zip(x, y):
        spl.AddPoint(i[0], i[1])
    return [spl.Evaluate(i) for i in newx]


def spline_interp_k(x, y, newx):
    spl = vtk.vtkKochanekSpline()
    spl.SetLeftConstraint(2)
    spl.SetRightConstraint(2)
    spl.SetDefaultTension(0.0)
    spl.SetDefaultContinuity(0.2)
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
        spline = vtk.vtkParametricSpline()
        spline.SetPoints(self.poly.GetPoints())
        spline.SetLeftConstraint(2)
        spline.SetRightConstraint(2)
        p, du = [0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]
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
            for j in zip(
                rel_dist_from_te,
                [
                    1 if k < len(rel_dist_from_te) // 2 else -1
                    for k in range(len(rel_dist_from_te))
                ],
            ):
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
        if len(self.points) < 3:
            return 0.0, 0.0
        g = list(zip(*self.points))
        return np.mean(g[1]), np.mean(g[2])

    def splits(self, r, r_relative):
        if not self.points:
            return (0, 0)
        out = (self.splines[0].Evaluate(r), self.splines[1].Evaluate(r))
        self.evaluations[int(round(r * 1e2) * 10)] = [out]
        return out

    def _find_top_and_bottom_points(self, webmesh):
        if not self.points:
            return
        rad = webmesh.GetPointData().GetArray("radius")
        rel_dist = webmesh.GetPointData().GetArray("d_rel_dist_from_te")
        for i in range(webmesh.GetNumberOfPoints()):
            rm = rad.GetValue(i)
            if self.web_root <= rm <= self.web_tip:
                rd = rel_dist.GetValue(i)
                pnt = webmesh.GetPoint(i)
                rmm = int(round(rm * 1e2) * 10)
                if rmm in self.evaluations and (
                    abs(rd - self.evaluations[rmm][0][0]) < 1e-6
                    or abs(rd - self.evaluations[rmm][0][1]) < 1e-6
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
        if not hasattr(self, "webmesh"):
            return
        self.webmesh.save(vtpfile)
        logger.info(f"Wrote mesh to {vtpfile}")

    def mesh(self, webmesh, n_cells):
        self._find_top_and_bottom_points(webmesh)
        points, pdata = self._create_points(n_cells)
        logger.info(f"Web {self.name}: created {len(points)} points")
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
        workdir = Path(filename).parent
        for i in self.webs:
            if hasattr(i, "webmesh"):
                i.write_mesh(str(workdir / f"{i.name}.vtp"))
