import vtk
import numpy as np
import math
import os
import pyvista as pv


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
