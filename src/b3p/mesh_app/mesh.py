#! /usr/bin/env python3
"""Mesh building functions for mesh_app."""

import os
import logging
import pickle
from copy import deepcopy as dc
import numpy as np
from pathlib import Path
import pyvista as pv
import vtk

logger = logging.getLogger(__name__)

# Copied from blade_section.py
class section:
    def __init__(self, x, y):
        self.x, self.y = x, y
        pnts = np.column_stack((x, y, np.zeros_like(x)))
        cells = np.column_stack(
            (
                2 * np.ones(len(x)).astype(int),
                np.arange(0, len(x)),
                np.arange(1, (len(x) + 1)) % len(x),
            )
        ).flatten()
        self.polydata = pv.PolyData(pnts, lines=cells)

    def local_to_global(self):
        """Transform the section to global coordinates from airfoil coordinates"""
        self.polydata.points = np.array(
            [self.polydata.points[:, i] for i in [1, 0, 2]]
        ).T

    def get_max_thickness(self, web_angle=0, n_points=50):
        """
        get the location (x) where the airfoil is thickest, used to offset the
        section for maximum building height.
        """
        bounds = self.polydata.GetBounds()
        dx = bounds[1] - bounds[0]
        px = np.linspace(bounds[0] + 0.15 * dx, bounds[1] - 0.4 * dx, n_points)
        plane = vtk.vtkPlane()
        clip = vtk.vtkCutter()

        clip.SetInputData(self.polydata)
        plane.SetNormal(
            np.cos(math.radians(web_angle)), math.sin(math.radians(web_angle)), 0
        )
        t = []
        tb = []
        for i in px:
            plane.SetOrigin(i, 0.5 * (bounds[2] + bounds[3]), 0)
            clip.SetCutFunction(plane)
            clip.Update()
            section = clip.GetOutput()
            top, bot = section.GetPoint(0), section.GetPoint(1)
            t.append(section.GetPoint(0)[1] - section.GetPoint(1)[1])
            tb.append((top, bot))

        return px[t.index(max(t))]

    def scale(self, scalefactor):
        self.polydata.scale(scalefactor, inplace=True)

    def twist(self, rz):
        self.polydata.rotate_z(rz, inplace=True)

    def translate(self, dx, dy, dz):
        self.polydata.translate([dx, dy, dz], inplace=True)

    def get_point(self, xy):
        return self.polydata.GetPoint(self.polydata.FindPoint((xy[0], xy[1], 0.0)))

    def get_pointlist(self, z_rotation=0):
        output = self.polydata.rotate_z(z_rotation, inplace=False)
        return output.points

    def to_xfoil(self, fname):
        if not os.path.isdir(os.path.dirname(fname)):
            os.makedirs(os.path.dirname(fname))
        with open(fname, "w") as f:
            for i in zip(self.x, self.y):
                f.write("%f      %f\n" % (i[0], i[1]))

    def get_te(self):
        te1 = self.polydata.GetPoint(0)
        te2 = self.polydata.GetPoint(self.polydata.GetNumberOfPoints() - 1)
        return (
            [0.5 * (i[0] + i[1]) for i in zip(te1, te2)],
            vtk.vtkMath.Distance2BetweenPoints(te1, te2) ** 0.5,
        )

# Copied from loft_utils.py
def load(fl, normalise=False):
    """
    load airfoil

    args:
        fl (str): filename

    kwargs:
        normalise (bool): flag determining whether the airfoil is normalised to
        unit length

    returns:
        [(x,y),...] airfoil point coordinates

    """
    d = []
    with open(fl, "r") as f:
        for i in f:
            with contextlib.suppress(Exception):
                xy = [float(j) for j in i.split()]
                if len(xy) in {2, 3}:
                    d.append(xy)
    x, y = list(zip(*d))
    if normalise:
        mx = min(x)
        dx = max(x) - min(x)
        x = [(i - mx) / dx for i in x]
        y = [i / dx for i in y]

    return list(zip(x, y))

# Copied from loft_utils.py
def interp(x, points):
    """
    Interpolate airfoil points using 3D vtkParametricSpline

    args:
        x : List of points in range [0,1]
        points: List of 3d points [(x,y,z),...] to interpolate through
    """
    pnts = vtk.vtkPoints()
    for i in points:
        pnts.InsertNextPoint(i[0], i[1], 0 if len(i) == 2 else i[2])
    spline = vtk.vtkParametricSpline()
    spline.SetPoints(pnts)
    spline.SetLeftConstraint(3)
    spline.SetLeftValue(1.0)
    spline.SetRightConstraint(3)
    spline.SetRightValue(1.0)
    p, du = [0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]
    spline.DerivativesAvailableOn()
    out = []
    for i in x:
        p[0] = i
        u = [0, 0, 0]
        spline.Evaluate(p, u, du)
        out.append(u)
    return list(zip(*out))

# Copied from splining.py
def intp_c(x, points, const=2, clamp=True):
    try:
        sc = vtk.vtkCardinalSpline()
        sc.SetLeftConstraint(const)
        sc.SetRightConstraint(const)
        if clamp:
            sc.ClampValueOn()
        else:
            sc.ClampValueOff()
        for i in points:
            sc.AddPoint(i[0], i[1])

        o = [sc.Evaluate(i) for i in x]
        return x, o

    except Exception as e:
        logger.error(f"Error occurred while interpolating points: {points}")
        logger.exception(e)

class blade:
    def __init__(
        self,
        chord,
        thickness,
        twist,
        dx,
        dy,
        z,
        airfoils,
        chordwise_sampling,
        np_spanwise=100,
    ):
        self.np_spanwise = np_spanwise
        self.np_chordwise = len(chordwise_sampling)
        self._load_airfoils(airfoils, chordwise_sampling)
        self._interpolate_planform(chord, thickness, twist, dx, dy, z)
        self._place_airfoils()

    def to_table(self, x, prefix="prebend_out"):
        cols = [
            "relative_r",
            "z",
            "prebend",
            "chord",
            "relative_thickness",
            "absolute_thickness",
            "twist",
        ]
        df = pd.DataFrame(
            np.array(
                [x]
                + [
                    np.interp(x, i[0], i[1])
                    for i in [
                        self.z,
                        self.dx,
                        self.chord,
                        self.thickness,
                        self.absolute_thickness,
                        self.twist,
                    ]
                ]
            ).T,
            columns=cols,
        )
        df.to_csv(f"{prefix}.csv", index=False, sep=";")

    def _load_airfoils(self, airfoils, x):
        self.airfoils = {}
        for i in sorted(airfoils):
            if isinstance(airfoils[i], str):
                if airfoils[i].find("du") != -1:
                    t = load(airfoils[i], normalise=True)
                else:
                    t = load(airfoils[i], normalise=False)
            else:
                t = airfoils[i]["xy"]
            self.airfoils[i] = interp(x, t)[:2]

    def _interpolate_planform(self, chord, thickness, twist, dx, dy, z):
        self.x = np.linspace(0, 1.0, self.np_spanwise)
        (
            self.input_chord,
            self.input_thickness,
            self.input_twist,
            self.input_dx,
            self.input_dy,
        ) = (
            list(zip(*dc(chord))),
            list(zip(*dc(thickness))),
            list(zip(*dc(twist))),
            list(zip(*dc(dx))),
            list(zip(*dc(dy))),
        )

        self.chord = self._intp_c(self.x, chord)
        self.twist = self._intp_c(self.x, twist)
        self.thickness = self._intp_c(self.x, thickness)
        self.dx = self._intp_c(self.x, dx)
        self.dy = self._intp_c(self.x, dy)
        self.z = self._intp_c(self.x, z)
        self.absolute_thickness = [
            self.x,
            [i[0] * i[1] for i in zip(self.chord[1], self.thickness[1])],
        ]

    def _intp_c(self, x, points):
        sc = vtk.vtkCardinalSpline()
        sc.SetLeftConstraint(2)
        sc.SetRightConstraint(2)
        sc.ClampValueOn()
        for i in points:
            sc.AddPoint(i[0], i[1])
        o = [sc.Evaluate(i) for i in x]
        return x, o

    def plot(self, fname="_out.png"):
        # Implement plotting if needed
        pass

    def _interpolate_airfoils(self):
        v = [
            list(
                zip(
                    [i for _ in self.airfoils[i][0]],
                    self.airfoils[i][0],
                    self.airfoils[i][1],
                )
            )
            for i in sorted(self.airfoils)
        ]
        nv = []
        for i in zip(*v):
            t, x, y = zip(*i)
            nx = np.interp(self.thickness[1], t, x)
            ny = np.interp(self.thickness[1], t, y)
            nv.append(list(zip(nx, ny)))

        sections = []
        for i in zip(*nv):
            x, y = list(zip(*i))
            sections.append(section(x, y))

        return sections

    def _place_airfoils(self):
        self.sections = self._interpolate_airfoils()

        twist_center = 0.5

        for i in zip(self.x, self.chord[1], self.twist[1], self.sections):
            i[3].translate(-twist_center, 0.0, 0.0)
            i[3].twist(i[2])
            i[3].scale((i[1], i[1], 1.0))

        for i in self.sections:
            i.local_to_global()

        for i in zip(self.sections, self.dx[1], self.dy[1], self.z[1]):
            i[0].translate(i[1], i[2], i[3])

    def export_variables(self, fname):
        var = {
            "dx": self.dx,
            "dy": self.dy,
            "z": self.z,
            "twist": self.twist,
            "chord": self.chord,
            "thickness": self.thickness,
            "absolute_thickness": self.absolute_thickness,
        }
        vv = {}
        for i in var:
            vv[i] = np.array(var[i]).tolist()

        with open(fname, "w") as f:
            json.dump(vv, f)

    def dump(self, fname="__sections.txt", z_rotation=0.0):
        fname = str(fname)  # Convert Path to string
        lst = [i.get_pointlist(z_rotation=z_rotation) for i in self.sections]
        if fname.endswith(".txt"):
            with open(fname, "wb") as f:
                f.write(str(lst).encode("utf-8"))
        elif fname.endswith(".pck"):
            with open(fname, "wb") as f:
                pickle.dump(lst, f)

    def export_xfoil(self, prefix="airfoil_out/_xf"):
        # Implement as needed
        pass

    def mesh(self, fname=None):
        """Generate the 3D mesh from the blade sections."""
        n_points = self.np_chordwise
        points = []
        for i in self.sections:
            for j in range(i.polydata.GetNumberOfPoints()):
                pt = i.polydata.GetPoint(j)
                points.append([pt[0], pt[1], pt[2]])

        points = np.array(points)

        cells = []
        for i in range(1, len(self.sections)):
            s0 = range((i - 1) * n_points, i * n_points)
            s1 = range(i * n_points, (i + 1) * n_points)
            cells.extend(
                [4, s0[j], s1[j], s1[(j + 1) % n_points], s0[(j + 1) % n_points]]
                for j in range(n_points)
            )

        self.poly = pv.PolyData(points, np.hstack(cells))
        if fname is not None:
            logger.info(f"Saving mesh to {fname}")
            self.poly.save(fname)


def optspace(n_points, base=0.2):
    """Alternative to linspace for sampling."""
    x = np.linspace(0, 4.0 * np.pi, n_points)
    sp = 1.0 + base - np.cos(x)
    x1 = np.array([sum(sp[:i]) for i in range(len(x))])
    x1 = x1 / max(x1)
    return x1

# Copied from webs.py
def write_web(
    loc,
    normal,
    mesh,
    name,
    rootcut=0.0,
    tipcut=100.0,
    tip=0.0,
    zone="d_rel_dist_from_te",
):
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
    if c is None or p is None or cc is None:
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
        return []
    lw = list(zip(*lw))
    ww = list(zip(*ww))
    rt = list(zip(*rta))
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
    return out


def build_webs(mesh, webs, prefix="__dum"):
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
        )
        web_meshes[i] = fea_web
    return web_meshes

# Copied from mesh_from_loft.py
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
        web(
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
        sec = geometry_section(r, r_rel, i, open_te=False)
        nsec.append(sec)
    blade = blade_shape(
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


def build_blade_mesh(config, workdir):
    pln = config["planform"]
    radii = np.linspace(0, 100, 100)
    web_inputs = config["mesh"]["webs"]
    base_vtp = workdir / "blade_geometry.vtp"
    web_intersections = build_webs(str(base_vtp), web_inputs, prefix="blade")
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
