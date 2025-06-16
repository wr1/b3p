import logging
import pandas as pd
import numpy as np
from copy import deepcopy as dc
from matplotlib import pyplot as plt
import pickle
import pyvista as pv
import json
from b3p.geometry import splining
from b3p.geometry import loft_utils
from b3p.geometry import blade_section

# Configure logging
logger = logging.getLogger(__name__)


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
        logger.info("Loading airfoils")
        self.airfoils = {}
        for i in sorted(airfoils):
            if type(airfoils[i]) == str:
                if airfoils[i].find("du") != -1:
                    logger.debug(f"Loading {airfoils[i]} normalized")
                    t = loft_utils.load(airfoils[i], normalise=True)
                else:
                    logger.debug(f"Loading {airfoils[i]} unnormalized")
                    t = loft_utils.load(airfoils[i], normalise=False)
            else:
                logger.debug(f"Loading airfoil {airfoils[i]['name']} at thickness {i}")
                t = airfoils[i]["xy"]

            self.airfoils[i] = loft_utils.interp(x, t)[:2]

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

        self.chord = splining.intp_c(self.x, chord)
        self.twist = splining.intp_c(self.x, twist)
        self.thickness = splining.intp_c(self.x, thickness)
        self.dx = splining.intp_c(self.x, dx)
        self.dy = splining.intp_c(self.x, dy)
        self.z = splining.intp_c(self.x, z)
        self.absolute_thickness = [
            self.x,
            [i[0] * i[1] for i in zip(self.chord[1], self.thickness[1])],
        ]

    def plot(self, name="_dum", fname="_out.png"):
        fig, ax = plt.subplots(3, 3, figsize=(15, 15))

        ax[0, 0].plot(self.x, self.chord[1], label=name)
        ax[0, 0].plot(
            self.input_chord[0], self.input_chord[1], "o", label=f"{name}_input"
        )
        ax[0, 0].set_xlabel("rel span (-)")
        ax[0, 0].legend(loc="best").get_frame().set_alpha(0.5)

        ax[2, 1].plot(self.x, self.chord[1], label=name)
        ax[2, 1].plot(
            self.input_chord[0], self.input_chord[1], "o", label=f"{name}_input"
        )
        ax[2, 1].grid(True)
        ax[2, 1].set_xlim(0.9, 1)

        ax[0, 1].plot(self.x, self.twist[1], label=name)
        ax[0, 1].plot(
            self.input_twist[0], self.input_twist[1], "o", label=f"{name}_input"
        )
        ax[0, 2].plot(self.x, self.thickness[1], label=name)
        ax[0, 2].plot(
            self.input_thickness[0],
            self.input_thickness[1],
            "o",
            label=f"{name}_input",
        )
        ax[1, 0].plot(
            self.absolute_thickness[0], self.absolute_thickness[1], label=name
        )
        ax[2, 0].plot(
            self.absolute_thickness[0], self.absolute_thickness[1], label=name
        )
        ax[2, 0].set_xlim(0, 0.2)

        ax[1, 1].plot(self.dx[0], self.dx[1], label=f"{name}_x")
        ax[1, 1].plot(self.input_dx[0], self.input_dx[1], "o", label=f"{name}_input")
        ax[1, 2].plot(self.dy[0], self.dy[1], label=f"{name}_y")
        ax[1, 2].plot(self.input_dy[0], self.input_dy[1], "o", label=f"{name}_input")

        titles = [
            "chord rotor_diam=%.3f" % (2.0 * max(self.z[1])),
            "twist",
            "relative thickness",
            "absolute thickness",
            "dx",
            "dy",
            "abs thickness root",
            "tip chord",
            "",
        ]

        for i in zip(ax.flatten(), titles):
            i[0].grid(True)
            i[0].set_title(i[1])

        ax[2, 2].remove()

        fig.tight_layout()
        fig.savefig(fname, dpi=100)

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
            sections.append(blade_section.section(x, y))

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

        json.dump(vv, open(fname, "w"))
        logger.info(f"Saved variables to {fname}")
        return var

    def dump(self, fname="__sections.txt", z_rotation=0.0):
        lst = [i.get_pointlist(z_rotation=z_rotation) for i in self.sections]
        if fname.endswith(".txt"):
            open(fname, "wb").write(str(lst).encode("utf-8"))
        elif fname.endswith(".pck"):
            pickle.dump(lst, open(fname, "wb"))

    def export_xfoil(self, prefix="airfoil_out/_xf"):
        for i in zip(self.sections, self.thickness[1], self.z[1]):
            nm = prefix + "_t_%.3f_r_%.3f" % (i[1], i[2])
            i[0].to_xfoil(nm.replace(".", "_") + ".dat")

    def mesh(self, fname=None):
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
