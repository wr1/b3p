from b3p import geom_utils
import vtk
import numpy as np
import math
import pyvista as pv

# from matplotlib import pyplot as plt
import copy


class section:
    """
    represents the geometry of a blade section

    args:
        r (list): list of radius position

        r_relative (list): list of relative radius positions

        points (list): list of points for this section

        min_te_thickness (float):  min TE thickness

    """

    def __init__(self, r, r_relative, points, min_te_thickness=0.002, open_te=True):
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
        """
        Open the trailing edge

        parameters:
            points (list): points
            te_thickness (float): trailing edge thickness
        """

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
        tol = 1e-9
        if webs == []:
            return np.linspace(0.0, 1.0, n_points)
        # compute where the web splits are on average (assuming they don't
        # cross)
        avspl = []
        for i in webs:
            avspl.extend(i.average_splits())

        # compute spacing between the splits
        intervals = sorted([0, 1] + list(avspl))
        isize = [intervals[i + 1] - intervals[i] for i in range(len(intervals) - 1)]

        for ps in panel_mesh_scale:
            if ps[0] < len(isize):
                isize[ps[0]] *= ps[1]

        # compute the number of cells there should be between each web split
        inum = [int(round(n_points * i / sum(isize))) for i in isize]

        if mult:
            isize[0] *= 2
            isize[-1] *= 2

        for _ in range(5):
            # consistent manner
            if (
                sum(inum) < n_points
            ):  # if not enough points, add one to the section with fewest points
                inum[inum.index(min(inum))] += 1
            elif (
                sum(inum) > n_points
            ):  # if too many, remove one from the section with most points
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
        """
        re-evaluate the spline using n_points points, when there are webs they
        are put as hard points into the point list
        Add a whole bunch of related coordinates to the points, that can later
        be used to drape plies
        """
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
                last = copy.deepcopy(u)
            distance_along_airfoil += geom_utils.distance(u, last)
            dist_from_te.append(geom_utils.distance(u, te))
            rel_dist_from_te.append(i)
            dist.append(distance_along_airfoil)
            x_abs.append(u[0])
            y_abs.append(u[1])
            last = copy.deepcopy(u)

        # compute a coordinate from the thickest point
        x, y, z = zip(*out)
        dist_miny = dist[y.index(min(y))]

        dist_sparcap = [i - dist_miny for i in dist]

        max_dist_from_te = max(dist_from_te)

        le_point = dist_from_te.index(max_dist_from_te)

        # chordline vector
        cline = [i[1] - i[0] for i in zip(out[le_point], te)]

        nor = [0, 0, 0]
        # get the in-plane normal to the chordline
        vtk.vtkMath.Cross(cline, [0, 0, -1], nor)
        chord = []
        for i in out:
            v1 = cline  # chord line
            v2 = nor  # normal to chord line
            v3 = [j[1] - j[0] for j in zip(te, i)]  # vector from TE to point

            chord.append(
                max(
                    0,
                    -(v3[0] / v1[0] - (v2[0] * v3[1]) / (v2[1] * v1[0]))
                    / (1.0 - (v2[0] * v1[1]) / (v2[1] * v1[0])),
                )
            )

        lechord = [-i + 1.0 for i in chord]
        lechord_absolute = [i * (cline[0] ** 2 + cline[1] ** 2) ** 0.5 for i in lechord]

        le_datum = [math.fabs(dist[le_point] - i) for i in dist]

        te_datum = [
            (dist[i] if i < le_point else max(dist) - dist[i]) for i in range(len(dist))
        ]

        # coordinates specific to the leading edge mould split line, defined as
        # the point where the y-coordinate is minimal
        sl_point = y_abs.index(min(y_abs))
        sl_dist = np.array(dist) - dist[sl_point]

        # suction and pressure sides have different false conditions to
        # distinguish them in case of centroid datasets
        suction_side = [1 if i < sl_point else -1 for i in range(len(dist))]
        pressure_side = [0 if i < sl_point else 1 for i in range(len(dist))]

        # add a coordinate that specifies the distance to the web attachment,
        # this is convenient for the specification of girders
        web_datums = {}
        for i in webs:
            splits = sorted(i.splits(self.r, self.r_relative))
            datum, datum_r = [], []
            for j in zip(rel_dist_from_te, pressure_side):
                datum.append(
                    (j[0] - splits[j[1]]) * distance_along_airfoil * (-1 if j[1] else 1)
                )
                datum_r.append((j[0] - splits[j[1]]) * (-1 if j[1] else 1))

            web_datums[f"d_{i.coordinate}"] = datum
            web_datums[f"d_{i.coordinate}_r"] = datum_r

        r = [self.r for _ in dist]
        # d_along_airfoil: distance from the TE to the point at hand along the
        # section
        # radius: radius position of the point
        # d_from_te: distance between current point and TE in a straight line
        # rel_dist_from_te: distance of point from TE scaled from 0 to 1
        # chord_length: chord length of the section the point is a part of
        # loc_le: d_along_airfoil coordinate of the leading edge point

        mdist = max(dist)
        datums = {
            "d_te": te_datum,
            "radius": r,
            "d_rel_dist_from_te": rel_dist_from_te,
            "d_abs_dist_from_te": dist,
            "d_abs_dist_from_bte": [-i + mdist for i in dist],
            "chord_length": [max_dist_from_te for _ in dist_from_te],
            "d_miny": dist_sparcap,
            "d_le": le_datum,
            "zone_ss": suction_side,
            "zone_ps": pressure_side,
            "d_te_r": [i / max(te_datum) for i in te_datum],
            "d_le_r": [i / max(le_datum) for i in le_datum],
            "d_chord": [1.0 - i for i in chord],
            "d_lechord": lechord,
            "d_lechord_abs": lechord_absolute,
            "d_x": x_abs,
            "d_y": y_abs,
            "d_sl": sl_dist,
            "d_sla": abs(sl_dist),
            "is_web": [0 for _ in r],
        }
        for i in web_datums:
            datums[i] = np.array(web_datums[i]).astype(np.float32)

        for i in added_datums.items():
            # print(self.r_relative, i[1][1], i[1][2])
            offs = np.interp(self.r_relative, i[1][1], i[1][2])
            datums[i[0]] = np.array(np.array(datums[i[1][0]]) + offs).astype(np.float32)

        return out, datums
