#! /usr/bin/env python3
"""Geometry building functions for geom_app."""

import os
from copy import deepcopy as dc
import numpy as np
import pandas as pd
import json
from pathlib import Path

# Inline or copy necessary functions from b3p.geometry
# Assuming we copy the blade class and related functions here to avoid imports

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
        # Simplified, assuming airfoils are provided
        self.airfoils = airfoils
        # Implement as needed

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
        # Simplified interpolation
        return [x, np.interp(x, [p[0] for p in points], [p[1] for p in points])]

    def plot(self, fname="_out.png"):
        # Implement plotting if needed
        pass

    def _interpolate_airfoils(self):
        # Simplified
        return []

    def _place_airfoils(self):
        self.sections = self._interpolate_airfoils()

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
        # Implement as needed
        pass

    def export_xfoil(self, prefix="airfoil_out/_xf"):
        # Implement as needed
        pass

    def mesh(self, fname=None):
        # Simplified mesh generation
        # Assume pyvista is available
        import pyvista as pv
        points = np.random.rand(100, 3)  # Placeholder
        cells = np.arange(100)  # Placeholder
        self.poly = pv.PolyData(points, cells)
        if fname:
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
    blade_obj.mesh(workdir / "blade_geometry.vtu")
    blade_obj.dump(workdir / "blade_geometry.pck")
    blade_obj.export_variables(workdir / "blade_geometry_variables.json")
    blade_obj.export_xfoil(
        prefix=os.path.join(workdir, "airfoil_out"),
    )
    blade_obj.plot(fname=workdir / "blade_geometry.png")
    n_sections = 50
    blade_obj.to_table(np.linspace(0, 1, n_sections), workdir / "blade_geometry_sca_50.csv")

    return blade_obj
