#! /usr/bin/env python3
"""Geometry building functions for geom_app."""

import os
import logging
import contextlib
import math
import pickle
from copy import deepcopy as dc
import numpy as np
import pandas as pd
import json
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


def build_blade_geometry(config, workdir):
    """Perform blade 3d model generation based on b3p dictionary."""
    pln = config["planform"]
    blade_obj = blade(
        pln["chord"],
        pln["thickness"],
        pln["twist"],
        dc(pln["dx"]),
        pln["dy"],
        pln["z"],
        config["aero"]["airfoils"],
        chordwise_sampling=optspace(
            config["planform"]["npchord"]
        ),
        np_spanwise=config["planform"]["npspan"],
    )

    # Use fixed file names in workdir
    blade_obj.mesh(workdir / "blade_geometry.vtp")
    blade_obj.dump(workdir / "blade_geometry.pck")
    blade_obj.export_variables(workdir / "blade_geometry_variables.json")
    blade_obj.export_xfoil(
        prefix=os.path.join(workdir, "airfoil_out"),
    )
    blade_obj.plot(fname=workdir / "blade_geometry.png")
    n_sections = 50
    blade_obj.to_table(np.linspace(0, 1, n_sections), workdir / "blade_geometry_sca_50.csv")

    return blade_obj
