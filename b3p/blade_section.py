import vtk
import numpy as np
import math
import os


class section:
    def __init__(self, x, y):
        self.x, self.y = x, y
        pnts = vtk.vtkPoints()
        for i in zip(x, y):
            pnts.InsertNextPoint(i[0], i[1], 0.0)
        self.polydata = vtk.vtkPolyData()
        self.polydata.SetPoints(pnts)

        cells = vtk.vtkCellArray()
        for i in range(1, pnts.GetNumberOfPoints()):
            cells.InsertNextCell(2)
            cells.InsertCellPoint(i - 1)
            cells.InsertCellPoint(i)
        self.polydata.SetLines(cells)

    def local_to_global(self):
        pnts = vtk.vtkPoints()
        for i in range(self.polydata.GetNumberOfPoints()):
            pnt = self.polydata.GetPoint(i)
            new_point = (pnt[1], pnt[0], pnt[2])
            pnts.InsertNextPoint(new_point)

        self.polydata.SetPoints(pnts)

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
        transform = vtk.vtkTransform()
        transform.Scale(scalefactor)
        transformfilter = vtk.vtkTransformFilter()
        transformfilter.SetTransform(transform)
        transformfilter.SetInputData(self.polydata)
        transformfilter.Update()
        self.polydata = transformfilter.GetOutput()

    def twist(self, rz):
        transform = vtk.vtkTransform()
        transform.RotateZ(rz)
        transformfilter = vtk.vtkTransformFilter()
        transformfilter.SetTransform(transform)
        transformfilter.SetInputData(self.polydata)
        transformfilter.Update()
        self.polydata = transformfilter.GetOutput()

    def translate(self, dx, dy, dz):
        transform = vtk.vtkTransform()
        transform.Translate(dx, dy, dz)
        transformfilter = vtk.vtkTransformFilter()
        transformfilter.SetTransform(transform)
        transformfilter.SetInputData(self.polydata)
        transformfilter.Update()
        self.polydata = transformfilter.GetOutput()

    def get_point(self, xy):
        return self.polydata.GetPoint(self.polydata.FindPoint((xy[0], xy[1], 0.0)))

    def get_pointlist(self, z_rotation=0):
        transform = vtk.vtkTransform()
        transform.RotateZ(z_rotation)
        transformfilter = vtk.vtkTransformFilter()
        transformfilter.SetTransform(transform)
        transformfilter.SetInputData(self.polydata)
        transformfilter.Update()
        output = transformfilter.GetOutput()
        return [output.GetPoint(i) for i in range(output.GetNumberOfPoints())]

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
