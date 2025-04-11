#! /usr/bin/env python

import vtk
import numpy as np
import json


def write_web(
    loc,
    normal,
    mesh,
    name,
    rootcut=0.0,
    tipcut=100.0,
    tip=0.0,
    zone="d_rel_dist_from_te",
):
    """Slices a plane through a mesh, and reports back a relative
    coordinate of top and bottom lines. This is used to represent geometrically straight entities
    in 3D as coordinates in local systems defined per section.

    Args:
        loc (np.ndarray): location of the plane
        normal (tuple): normal of the plane
        mesh (str): mesh to slice
        name (str): name of the web
        rootcut (float, optional): cut the root of the web. Defaults to 0.0.
        tipcut (float, optional): cut the tip of the web. Defaults to 100.0.
        tip (float, optional): tip of the blade. Defaults to 0.0.
        zone (str, optional): zone of the blade. Defaults to "d_rel_dist_from_te".

    Returns:
        list: A list of tuples, each containing (rr, wwl, lwl, rt1[i])
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

        # Prepare data for JSON serialization
        data = {
            "name": name,
            "data": out,
            "points": {"lwp": [list(p) for p in lwp], "wwp": [list(p) for p in wwp]},
        }

        # Write data to a JSON file
        with open(f"{name}.json", "w") as f:
            json.dump(data, f, indent=4)

    # open("%s.txt" % name, "wb").write(str(out).encode("utf-8"))
    # open("%s_points.txt" % name, "wb").write(str([lwp, wwp]).encode("utf-8"))

    return out


def build_webs(mesh, webs, prefix="__dum"):
    """Builds web meshes based on provided web definitions.
    Args:
        mesh (object): The base mesh object to which the webs will be attached.
        webs (dict): A dictionary defining the webs to be created.
            Each key in the dictionary represents the name of a web, and the
            corresponding value is a dictionary containing the web's properties,
            including:
                "origin" (list/array-like): The origin point of the web.
                "z_start" (float): The z-coordinate where the web starts.
                "z_follow_blade" (float): The z-coordinate where the web follows the blade.
                "z_end" (float): The z-coordinate where the web ends.
                "orientation" (list/array-like, optional): The normal vector
                    defining the orientation of the web. Defaults to (0, 1, 0).
        prefix (str, optional): A prefix to be added to the name of each web mesh.
            Defaults to "__dum".
    Returns:
        dict: A dictionary containing the generated web meshes. The keys of the
            dictionary are the names of the webs, and the values are the
            corresponding mesh objects.
    """

    web_meshes = {}
    for i in webs:
        normal = (0, 1, 0)
        name = str(prefix) + "_" + i

        if "orientation" in webs[i]:
            normal = webs[i]["orientation"]

        fea_web = write_web(
            np.array(webs[i]["origin"]),
            normal,
            mesh,
            name,
            rootcut=webs[i]["z_start"],
            tipcut=webs[i]["z_follow_blade"],
            tip=webs[i]["z_end"],
        )
        web_meshes[i] = fea_web

    return web_meshes
