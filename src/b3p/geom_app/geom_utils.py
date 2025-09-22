import vtk
import math


def spline_interp(x, y, newx):
    spl = vtk.vtkCardinalSpline()
    spl.SetLeftConstraint(2)
    spl.SetRightConstraint(2)
    for i in zip(x, y):
        spl.AddPoint(i[0], i[1])

    return [spl.Evaluate(i) for i in newx]


def spline_interp_k(x, y, newx):
    spl = vtk.vtkKochanekSpline()
    spl.SetLeftConstraint(2)
    spl.SetRightConstraint(2)
    spl.SetDefaultTension(0.0)
    spl.SetDefaultContinuity(0.2)
    for i in zip(x, y):
        spl.AddPoint(i[0], i[1])

    return [spl.Evaluate(i) for i in newx]


def distance(point1, point2):
    return math.sqrt(sum((i[1] - i[0]) ** 2 for i in zip(point1, point2)))
