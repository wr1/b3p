import logging
import vtk
import numpy as np
import pyvista as pv

logger = logging.getLogger(__name__)


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


class web:
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
        out = (0, 0)
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
        self.mesh.save(vtpfile)
        logger.info(f"Wrote mesh to {vtpfile}")

    def mesh(self, mesh, n_cells):
        self._find_top_and_bottom_points(mesh)
        points, pdata = self._create_points(n_cells)

        cells = self._create_quad_connectivity(n_cells, len(points), self.flip_normal)

        self.mesh = pv.PolyData(points, faces=cells)

        for i in pdata:
            self.mesh.point_data[i] = np.array(pdata[i]).astype(np.float32)
