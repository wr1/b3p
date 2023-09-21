from b3p import splining
from b3p import loft_utils
from b3p import blade_section

import pandas as pd
import numpy as np
from copy import deepcopy as dc
from matplotlib import pyplot as plt
import pickle
import math
import pyvista as pv


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
        """Build a blade geometry from a yaml file
        :param chord: list of chord lengths
        :param thickness: list of relative thicknesses
        :param twist: list of twist angles
        :param dx: list of x offsets
        :param dy: list of y offsets
        :param z: list of z offsets
        :param airfoils: dict of airfoils
        :param chordwise_sampling: list of chordwise sampling points
        :param np_spanwise: number of spanwise points"""
        self.np_spanwise = np_spanwise
        self.np_chordwise = len(chordwise_sampling)
        self._load_airfoils(airfoils, chordwise_sampling)
        self._interpolate_planform(chord, thickness, twist, dx, dy, z)
        self._place_airfoils()

    def to_table(self, x, prefix="prebend_out"):
        """
        write out a table of the blade

        parameters:
            x (list) : list of sections where airfoils are to be interpolated
            prefix (str) : prefix for output file
        """
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
        """
        fill self.airfoil by interpolating from input airfoil set

        args:
            airfoils (dict) : airfoils dict

            x (list) : list of sections where airfoils are to be interpolated
        """
        print("** loading airfoils")
        self.airfoils = {}
        for i in sorted(airfoils):
            if type(airfoils[i]) == str:
                if airfoils[i].find("du") != -1:
                    print(f"** load {airfoils[i]} normalised")
                    t = loft_utils.load(airfoils[i], normalise=True)
                else:
                    print(f"** load {airfoils[i]} unnormalised")
                    t = loft_utils.load(airfoils[i], normalise=False)  # fix so we don't
                    # normalise flatback root
            else:
                # load an airfoil from the self-contained format, which is a dict with keys xy
                print(f"** loading airfoil {airfoils[i]['name']} at thickness {i}")
                t = airfoils[i]["xy"]

            self.airfoils[i] = loft_utils.interp(x, t)[:2]

    def _interpolate_planform(self, chord, thickness, twist, dx, dy, z):
        """
        spline the planform based on parameters

        parameters:
            chord (list) : chord distribution [(r,chord),..]
            thickness (list) : thickness distribution [(r,t),..] --> relative
            twist (list) : twist distribution [(r,twist),..]
            dx (list) : x distribution [(r,x),..]
            dy (list) : y distribution [(r,y),..]
            z (list) : z distribution [(r,z),..]
            flatten_lw (bool): flag indicating whether to run LW side flattening
        """
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
        """
        plot parameter distributions

        args:
            name (str) : blade candidate name

            fname (str): plot output file name

        """
        plt.figure(figsize=(15, 15))
        plt.subplot(3, 3, 1)
        plt.title("chord rotor_diam=%.3f" % (2.0 * max(self.z[1])))
        plt.plot(self.x, self.chord[1], label=name)
        plt.plot(self.input_chord[0], self.input_chord[1], "o", label=f"{name}_input")
        plt.xlabel("rel span (-)")
        plt.ylabel("chord (m)")
        plt.grid(True)
        plt.legend(loc="best").get_frame().set_alpha(0.5)

        plt.subplot(3, 3, 8)
        plt.plot(self.x, self.chord[1], label=name)
        plt.plot(self.input_chord[0], self.input_chord[1], "o", label=f"{name}_input")
        plt.grid(True)
        plt.xlim(0.9, 1)

        plt.title("tip chord")
        plt.subplot(3, 3, 2)
        plt.plot(self.x, self.twist[1], label=name)
        plt.plot(self.input_twist[0], self.input_twist[1], "o", label=f"{name}_input")
        self.title_plot("twist", 3)
        plt.plot(self.x, self.thickness[1], label=name)
        plt.plot(
            self.input_thickness[0],
            self.input_thickness[1],
            "o",
            label=f"{name}_input",
        )
        self.title_plot("rel thickness", 4)
        plt.plot(self.absolute_thickness[0], self.absolute_thickness[1], label=name)
        self.title_plot("abs thickness", 7)
        plt.plot(self.absolute_thickness[0], self.absolute_thickness[1], label=name)
        plt.title("abs thickness")
        plt.grid(True)
        plt.xlim(0, 0.2)

        plt.subplot(3, 3, 5)
        plt.plot(self.dx[0], self.dx[1], label=f"{name}_x")
        plt.plot(self.input_dx[0], self.input_dx[1], "o", label=f"{name}_input")
        plt.plot(self.dy[0], self.dy[1], label=f"{name}_y")
        plt.plot(self.input_dy[0], self.input_dy[1], "o", label=f"{name}_input")
        plt.legend(loc="best").get_frame().set_alpha(0.5)
        plt.savefig(fname, dpi=100)

    # TODO Rename this here and in `plot`
    def title_plot(self, arg0, arg1):
        plt.title(arg0)
        plt.grid(True)
        plt.subplot(3, 3, arg1)

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

        # build the blade up out of sections in two loops
        # first, scale and twist the section
        # a variable is created, analogous to focus, which represents the
        # fraction of the chord around which the twist is defined
        twist_center = 0.5

        for i in zip(self.x, self.chord[1], self.twist[1], self.sections):
            i[3].translate(-twist_center, 0.0, 0.0)
            i[3].twist(i[2])
            i[3].scale((i[1], i[1], 1.0))

        for i in self.sections:
            i.local_to_global()

        # then, translate the section coordinates to global
        for i in zip(self.sections, self.dx[1], self.dy[1], self.z[1]):
            i[0].translate(i[1], i[2], i[3])

    def export_variables(self, fname):
        var = {
            "dx": self.dx,
            # "dxf": self.dxf,
            "dy": self.dy,
            "z": self.z,
            "twist": self.twist,
            "chord": self.chord,
            "thickness": self.thickness,
            "absolute_thickness": self.absolute_thickness,
        }
        open(fname, "w").write(str(var))
        return var

    def dump(self, fname="__sections.txt", z_rotation=0.0):
        "dump to a sections list for use in FEA (with webs)"
        lst = [i.get_pointlist(z_rotation=z_rotation) for i in self.sections]
        if fname.endswith(".txt"):
            open(fname, "wb").write(str(lst).encode("utf-8"))
        elif fname.endswith(".pck"):
            pickle.dump(lst, open(fname, "wb"))

    def export_xfoil(self, prefix="airfoil_out/_xf"):
        """Export the sections to xfoil format

        parameters
        ----------
        prefix : str
            prefix to use for the filenames
        """
        for i in zip(self.sections, self.thickness[1], self.z[1]):
            nm = prefix + "_t_%.3f_r_%.3f" % (i[1], i[2])
            i[0].to_xfoil(nm.replace(".", "_") + ".dat")

    def mesh(self, fname=None):
        """Join up the sections to make a polydata object and save it to a file

        parameters
        ----------
        fname : str
            filename to save the polydata to.  If None, don't save to file.
        """
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
        if fname != None:
            print("saving to", fname)
            self.poly.save(fname)
