from b3p import splining
from b3p import loft_utils
from b3p import blade_section


import numpy as np
from copy import deepcopy as dc
from matplotlib import pyplot as plt
import pickle
import math
import vtk
import json


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
        np_spanwise=100,
        chordwise_sampling=[],
        offset_optimal=True,
        interpolate_method=1,
        flatten_lw=True,
        barrel_length=0.0,
        offset_clamp_points=[0.32, 0.55, 0.7],
    ):
        """
        sequence of what needs to be done:

        #. load in airfoils

        #.  resample airfoils to fixed number of points around circumference
            using spline interpolation

        #. sample thickness, chord, twist and centerline at required z positions

        #. create linear interpolations of the airfoil table at sampled thickness
            values

        #. scale, rotate, offset sections to exist on the

        Args:
            chord : chord distribution [(r,chord),..]

            thickness : thickness distribution [(r,t),..] --> relative

            twist : twist distribution [(r,twist),..]

            centerline : centerline spline [(x,y,z),..]

            airfoils : table of airfoils (key is thickness) {.35:'naca35'}

            np_chordwise (int) : number of chordwise points

            np_spanwise (int): number of spanwise points

            chordwise_sampling (list) : chordwise list of t coordinates
                (overrides np_chordwise)

            offset_optimal (bool): flag indicating whether to run optimal
                offsetting of the sparcap

            interpolate_method (int): Type of interpolation (1 being linear, 2
            KochanekSpline , 3 CardinalSpline, 4 SmoothCurve )

            flatten_lw (bool): flag indicating whether to run LW side flattening

            barrel_length (float): length of root that stays cylindrical

            offset_clamp_points (list[float]): list of points closest to which
                the optimal offsetting points are taken

        """
        self.np_spanwise = np_spanwise
        self.barrel_length = barrel_length
        self.np_chordwise = len(chordwise_sampling)

        self._load_airfoils(airfoils, chordwise_sampling)
        self._interpolate_planform(
            chord, thickness, twist, dx, dy, z, flatten_lw=flatten_lw
        )
        self._interpolate_airfoils(
            offset_optimal=offset_optimal,
            interpolate_method=interpolate_method,
            offset_clamp_points=offset_clamp_points,
        )

    def to_table(self, prefix="prebend_out", x=[]):
        if x == []:
            x = self.dy[0]

        f = open("%s.csv" % prefix, "w")
        f.write(
            "relative_r;z;prebend;chord;relative_thickness;absolute_thickness;twist;dx\n"
        )

        z, dy, ch, th, thr, tw, dxf = [
            np.interp(x, i[0], i[1])
            for i in [
                self.z,
                self.dx,
                self.chord,
                self.thickness,
                self.absolute_thickness,
                self.twist,
                (
                    self.dy[0],
                    [i[1] - i[2] for i in zip(self.chord[1], self.dxf[1], self.dy[1])],
                ),
            ]
        ]

        for i in zip(x, z, dy, ch, th, thr, tw, dxf):
            f.write("%f;%f;%f;%f;%f;%f;%f;%f\n" % i)
        f.close()

    def _load_airfoils(self, airfoils, x=[]):
        """
        fill self.airfoil by interpolating from input airfoil set

        args:
            airfoils (dict) : airfoils dict

            x (list) : list of sections where airfoils are to be interpolated
        """
        self.airfoils = {}
        for i in sorted(airfoils):
            if type(airfoils[i]) == str:
                if airfoils[i].find("du") != -1:
                    print(f"load {airfoils[i]} normalised")
                    t = loft_utils.load(airfoils[i], normalise=True)
                else:
                    print(f"load {airfoils[i]} unnormalised")
                    t = loft_utils.load(airfoils[i], normalise=False)  # fix so we don't
                    # normalise flatback root
            else:
                # load an airfoil from the self-contained format, which is a dict with keys xy
                t = airfoils[i]["xy"]

            self.airfoils[i] = loft_utils.interp(x, t)[:2]

    def _interpolate_planform(
        self, chord, thickness, twist, dx, dy, z, flatten_lw=True
    ):
        """
        spline the planform based on parameters
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

        # here we take the thickness, and determine a y offset such that the
        # suction side sparcap runs flat for the first half of the blade
        dy_flat = [
            self.x,
            [
                0.5 * self.absolute_thickness[1][0] - 0.5 * i
                for i in self.absolute_thickness[1]
            ],
        ]

        if flatten_lw:
            mid_offset = dy_flat[1][int(len(dy_flat[1]) / 2)]
            dy_flat[1] = [
                i[1] - 2.0 * i[0] * mid_offset for i in zip(dy_flat[0], dy_flat[1])
            ]

            dynew = list(zip(dy_flat[0], dy_flat[1]))[: int(len(self.x) / 2)]

            for i in dy:
                if i[0] > 0.7:
                    dynew.append(i)

            if self.barrel_length > 0.00001:
                bpnts = []
                xm = 0.0
                for i in np.linspace(
                    0, self.barrel_length / (self.z[1][-1] - self.z[1][0]), 20
                ):
                    bpnts.append((i, 0))
                    xm = i

                for i in dynew:
                    if i[0] > 0.2:
                        bpnts.append(i)
                dynew = bpnts
            self.dy = splining.intp_c(self.x, dynew)

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
        self._extracted_from_plot_31("twist", 3)
        plt.plot(self.x, self.thickness[1], label=name)
        plt.plot(
            self.input_thickness[0],
            self.input_thickness[1],
            "o",
            label=f"{name}_input",
        )
        self._extracted_from_plot_31("rel thickness", 4)
        plt.plot(self.absolute_thickness[0], self.absolute_thickness[1], label=name)
        self._extracted_from_plot_31("abs thickness", 7)
        plt.plot(self.absolute_thickness[0], self.absolute_thickness[1], label=name)
        plt.title("abs thickness")
        plt.grid(True)
        plt.xlim(0, 0.2)

        plt.subplot(3, 3, 5)
        plt.plot(self.dx[0], self.dx[1], label="%s_x" % name)
        plt.plot(self.input_dx[0], self.input_dx[1], "o", label=f"{name}_input")
        plt.plot(self.dy[0], self.dy[1], label="%s_y" % name)
        plt.plot(self.input_dy[0], self.input_dy[1], "o", label=f"{name}_input")
        plt.legend(loc="best").get_frame().set_alpha(0.5)
        self._extracted_from_plot_31("xy offsets", 6)
        plt.plot(
            list(zip(*self.mx_thickness_loc))[0],
            list(zip(*self.mx_thickness_loc))[1],
            "o",
            label="control points",
        )
        plt.plot(self.x, self.mmxt, ".", label="sectionwise optimal")
        plt.plot(self.dxf[0], self.dxf[1], label="offset used")
        plt.ylabel("dist to 0.3 x chord (m)")
        plt.legend(loc="best")
        plt.grid(True)

        plt.savefig(fname, dpi=100)

    # TODO Rename this here and in `plot`
    def _extracted_from_plot_31(self, arg0, arg1):
        plt.title(arg0)
        plt.grid(True)

        plt.subplot(3, 3, arg1)

    def _interpolate_airfoils(
        self,
        interpolate_method=1,
        offset_optimal=True,
        offset_clamp_points=[0.32, 0.55, 0.7],
    ):
        v = []
        for i in sorted(self.airfoils):
            v.append(
                list(
                    zip(
                        [i for j in self.airfoils[i][0]],
                        self.airfoils[i][0],
                        self.airfoils[i][1],
                    )
                )
            )

        nv = []
        for i in zip(*v):
            t, x, y = zip(*i)
            if interpolate_method == 0:
                nx = np.interp(self.thickness[1], t, x)
                ny = np.interp(self.thickness[1], t, y)
            elif interpolate_method == 1:
                # interpolate along the length of the blade using a splining.intp_k
                dum, nx = splining.intp_k(self.thickness[1], zip(t, x))
                dum, ny = splining.intp_k(self.thickness[1], zip(t, y))
            elif interpolate_method == 2:
                dum, nx = splining.intp_c(self.thickness[1], zip(t, x))
                dum, ny = splining.intp_c(self.thickness[1], zip(t, y))
            elif interpolate_method == 3:
                dum, nx = splining.intp_sc(self.thickness[1], zip(t, x))
                dum, ny = splining.intp_sc(self.thickness[1], zip(t, y))
            else:
                exit("%i is not a valid interpolate_method" % i)
            nv.append(list(zip(nx, ny)))

        self.sections = []
        for i in zip(*nv):
            x, y = list(zip(*i))
            self.sections.append(blade_section.section(x, y))

        # build the blade up out of sections in two loops
        # first, scale and twist the section
        mx_thickness_loc = []
        c = 0

        # a variable is created, analogous to focus, which represents the
        # fraction of the chord around which the twist is defined
        twist_center = 0.3

        self.mmxt = []
        for i in zip(self.x, self.chord[1], self.twist[1], self.sections):
            # with -0.3 (since the section is not scaled, this is 0.3*chord)
            i[3].translate(-twist_center, 0.0, 0.0)
            i[3].twist(i[2])
            i[3].scale((i[1], i[1], 1.0))
            mxt = i[3].get_max_thickness()
            self.mmxt.append(mxt)
            for j in sorted(offset_clamp_points):
                if i[0] >= j:
                    mx_thickness_loc.append((i[0], mxt))
                    offset_clamp_points.remove(j)
                    break
            c += 1

        mx_thickness_loc.append((i[0], i[3].get_max_thickness()))

        dxf = list(splining.intp_c(self.x, mx_thickness_loc))

        dxf[1] = np.array(dxf[1])

        if offset_optimal == False:
            dxf[1] = np.zeros(len(dxf[1]))

        # use the local coordinate to offset the section so that the thickest
        # point lines up with the pitch axis
        if offset_optimal:
            fpx, fpy = [], []
            for i in zip(self.sections, dxf[1], self.twist[1]):
                fpx.append(math.sin(math.radians(i[2])) * i[1])
                fpy.append(-i[1] + (1.0 - np.cos(np.radians(i[2]))) * i[1])

        self.dx = (self.dx[0], np.array(self.dx[1]) - np.array(fpx))
        self.dy = (self.dy[0], np.array(self.dy[1]) + np.array(fpy))

        for i in self.sections:
            i.local_to_global()

        # then, translate the section coordinates to global
        for i in zip(self.sections, self.dx[1], self.dy[1], self.z[1]):
            i[0].translate(i[1], i[2], i[3])

        self.dxf = dxf
        self.mx_thickness_loc = mx_thickness_loc

    def export_variables(self, fname):
        var = {
            "dx": self.dx,
            "dxf": self.dxf,
            "dy": self.dy,
            "z": self.z,
            "twist": self.twist,
            "chord": self.chord,
            "thickness": self.thickness,
            "absolute_thickness": self.absolute_thickness,
        }
        json.dump(var, open(fname, "w"))
        # open(fname, "w").write(str(var))
        return var

    def dump(self, fname="__sections.txt", z_rotation=0.0):
        "dump to a sections list for use in FEA (with webs)"
        lst = []
        for i in self.sections:
            lst.append(i.get_pointlist(z_rotation=z_rotation))

        if fname.endswith(".txt"):
            open(fname, "wb").write(str(lst).encode("utf-8"))
        elif fname.endswith(".pck"):
            pickle.dump(lst, open(fname, "wb"))

    def export_xfoil(self, prefix="airfoil_out/_xf"):
        for i in zip(self.sections, self.thickness[1], self.z[1]):
            nm = prefix + "_t_%.3f_r_%.3f" % (i[1], i[2])
            i[0].to_xfoil(nm.replace(".", "_") + ".dat")

    def mesh(self, fname=""):
        "join up the sections"
        n_points = self.np_chordwise
        vp = vtk.vtkPoints()
        for i in self.sections:
            for j in range(i.polydata.GetNumberOfPoints()):
                pt = i.polydata.GetPoint(j)
                vp.InsertNextPoint((pt[0], pt[1], pt[2]))

        self.poly = vtk.vtkPolyData()
        self.poly.SetPoints(vp)

        cells = vtk.vtkCellArray()
        for i in range(1, len(self.sections)):
            s0 = range((i - 1) * n_points, i * n_points)
            s1 = range(i * n_points, (i + 1) * n_points)
            for j in range(n_points):
                cells.InsertNextCell(3)
                cells.InsertCellPoint(s0[j])
                cells.InsertCellPoint(s1[j])
                cells.InsertCellPoint(s1[(j + 1) % n_points])

                cells.InsertNextCell(3)
                cells.InsertCellPoint(s0[j])
                cells.InsertCellPoint(s1[(j + 1) % n_points])
                cells.InsertCellPoint(s0[(j + 1) % n_points])

        # bottom and top caps
        s0 = range(0, n_points)
        for i in range(0, int(math.ceil(0.5 * n_points) - 1)):
            t1 = (s0[i], s0[n_points - i - 2], s0[n_points - i - 1])
            t2 = (s0[i], s0[i + 1], s0[n_points - i - 2])

            for j in [t1, t2]:
                if j[1] != j[2]:
                    cells.InsertNextCell(3)
                    for k in j:
                        cells.InsertCellPoint(k)

        s0 = s1
        for i in range(0, int(math.ceil(0.5 * n_points) - 1)):
            t1 = (s0[i], s0[n_points - i - 1], s0[n_points - i - 2])
            t2 = (s0[i], s0[n_points - i - 2], s0[i + 1])

            for j in [t1, t2]:
                if j[1] != j[2]:
                    cells.InsertNextCell(3)
                    for k in j:
                        cells.InsertCellPoint(k)

        self.poly.SetPolys(cells)

        if fname != "":
            wr = vtk.vtkSTLWriter()
            wr.SetFileTypeToBinary()
            wr.SetFileName(fname)
            wr.SetInputData(self.poly)
            wr.Write()

            wr = vtk.vtkXMLPolyDataWriter()
            wr.SetFileName(fname.replace(".stl", ".vtp"))
            wr.SetInputData(self.poly)
            wr.Write()
