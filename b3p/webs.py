#! /usr/bin/env python

import vtk
import numpy as np


def write_web(
    loc,
    normal,
    mesh,
    rootcut=0.0,
    tipcut=100.0,
    tip=0.0,
    zone="d_rel_dist_from_te",
):
    """
    Slices a plane through a mesh, and reports back a relative
    coordinate of top and bottom lines. This is used to represent geometrically straight entities
    in 3D as coordinates in local systems defined per section

    :param loc: location of the plane
    :param normal: normal of the plane
    :param mesh: mesh to slice
    :param rootcut: cut the root of the web
    :param tipcut: cut the tip of the web
    :param tip: tip of the blade
    """

    rd = vtk.vtkXMLPolyDataReader()
    rd.SetFileName(mesh)
    rd.Update()
    poly = rd.GetOutput()

    plane = vtk.vtkPlane()
    plane.SetOrigin(loc)
    plane.SetNormal(normal)

    clip = vtk.vtkCutter()
    clip.SetCutFunction(plane)
    clip.SetInputData(poly)

    clip.Update()
    out = clip.GetOutput()
    points = []

    c = out.GetPointData().GetArray("d_rel_dist_from_te")
    cc = out.GetPointData().GetArray("d_abs_dist_from_te")
    p = out.GetPointData().GetArray(zone)

    lw, ww = [], []
    lwp, wwp = [], []
    rta = []
    minr, maxr = 1000, -1000
    for i in range(out.GetNumberOfPoints()):
        if c.GetValue(i) > 0.5:
            lw.append((out.GetPoint(i)[2], p.GetValue(i)))
            lwp.append(out.GetPoint(i))
            rta.append((out.GetPoint(i)[2], cc.GetValue(i) / c.GetValue(i)))
        else:
            ww.append((out.GetPoint(i)[2], p.GetValue(i)))
            wwp.append(out.GetPoint(i))
        minr = min(minr, out.GetPoint(i)[2])
        maxr = max(maxr, out.GetPoint(i)[2])
        points.append(out.GetPoint(i))

    lw = list(zip(*lw))
    ww = list(zip(*ww))
    rt = list(zip(*rta))
    # interpolate the results
    r = np.linspace(minr, maxr, 400)
    lw1 = np.interp(r, lw[0], lw[1])
    ww1 = np.interp(r, ww[0], ww[1])
    rt1 = np.interp(r, rt[0], rt[1])

    out = []
    for i, rr in enumerate(r):
        if rr <= tipcut:
            lwl, wwl = lw1[i], ww1[i]

        out.append((rr, wwl, lwl, rt1[i]))

        if i > 0 and r[i - 1] < rootcut and rr > rootcut:
            for j in range(i):
                out[j][1] = wwl
                out[j][2] = lwl

    out = sorted(out)
    if tip > out[-1][0]:
        out.append((tip, out[-1][1], out[-1][2], out[-1][3]))

    return out


def build_webs(mesh, webs, prefix="__dum"):
    web_meshes = {}
    for i in webs:
        normal = (0, 1, 0)

        if "orientation" in webs[i]:
            normal = webs[i]["orientation"]

        fea_web = write_web(
            np.array(webs[i]["origin"]),
            normal,
            mesh,
            rootcut=webs[i]["z_start"],
            tipcut=webs[i]["z_follow_blade"],
            tip=webs[i]["z_end"],
        )
        web_meshes[i] = fea_web

    return web_meshes
