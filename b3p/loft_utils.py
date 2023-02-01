#! /usr/bin/env python
import contextlib
import os
import numpy as np
import vtk


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
    print(f"loading airfoil {fl}")
    for i in open(fl, "r"):
        with contextlib.suppress(Exception):
            xy = [float(j) for j in i.split()]
            if len(xy) in {2, 3}:
                d.append(xy)
    x, y = list(zip(*d))
    if normalise:
        mx = min(x)
        dx = max(x) - min(x)
        print("normalise factor %f" % dx)
        x = [(i - mx) / dx for i in x]
        y = [i / dx for i in y]

    return list(zip(x, y))


def optspace(n_points, base=0.2):
    """
    alternative to linspace for sampling that puts more points near the TE and
    LE of an airfoil
    """
    lep = 0.5
    x = np.linspace(0, 4.0 * np.pi, n_points)
    sp = 1.0 + base - np.cos(x)
    x1 = np.array([sum(sp[:i]) for i in range(len(x))])
    x1 = x1 / max(x1)
    return x1


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
