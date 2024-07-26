#! /usr/bin/env python3

import sys
import vtk
import os

import multiprocessing
from functools import partial
import argparse
import numpy as np
import pandas as pd
from numpy import array
from scipy.spatial import distance
import math
import copy
import pyvista as pv


def get_eids(sec, n):
    g = vtk.vtkIdList()
    sec.GetPointCells(n, g)
    return [g.GetId(i) for i in range(g.GetNumberOfIds())]


def get_nids(sec, e):
    l = vtk.vtkIdList()
    sec.GetCellPoints(e, l)
    return [l.GetId(i) for i in range(l.GetNumberOfIds())]


def join_up(nodes, epid, stacks, sec, poly):
    """
    take nodes where ply drops happen, link them up using vtkDelaunay2D (using
    scaled down thicknesses
    """
    cells = poly.GetPolys()
    mat = poly.GetCellData().GetArray("mat")
    angle = poly.GetCellData().GetArray("angle")
    for inn in nodes:  # loop over all nodes that need joining up
        i, el = inn
        s0, s1 = np.array(stacks[el[0]]), np.array(stacks[el[1]])

        s0 /= np.linalg.norm(s0)
        s1 /= np.linalg.norm(s1)

        pnts = vtk.vtkPoints()
        for j in s0:
            pnts.InsertNextPoint(0, j, 0)
        for j in s1:
            pnts.InsertNextPoint(5, j, 0)

        pl = vtk.vtkPolyData()
        pl.SetPoints(pnts)
        d2 = vtk.vtkDelaunay2D()
        d2.SetInputData(pl)
        d2.Update()
        dd = d2.GetOutput()

        # new linking, fixes changes in ordering of elements
        x0 = 1
        ni0 = get_nids(sec, el[0])
        ni1 = get_nids(sec, el[1])
        if ni0[0] == ni1[1]:
            x0 = 0
        nmap = list(
            range(epid[el[0]] + x0, epid[el[0]] + x0 + (len(s0)) * 2, 2)
        ) + list(range(epid[el[1]] + 1 - x0, epid[el[1]] + 1 - x0 + (len(s1)) * 2, 2))
        for j in range(dd.GetNumberOfCells()):
            c = dd.GetCell(j)
            cells.InsertNextCell(4)
            for k in range(c.GetNumberOfPoints()):
                cells.InsertCellPoint(nmap[c.GetPointId(k)])
            cells.InsertCellPoint(nmap[c.GetPointId(k)])
            mat.InsertNextTuple1(-1)
            angle.InsertNextTuple1(0)


def align_normals(align):
    for i in range(align.GetNumberOfCells()):
        c = align.GetCell(i)
        p = c.GetPoints()
        if p.GetNumberOfPoints() > 2:
            normal = [0, 0, 0]
            vtk.vtkTriangle().ComputeNormal(
                p.GetPoint(0), p.GetPoint(1), p.GetPoint(2), normal
            )
            if normal[2] > 0:
                align.ReverseCell(i)


def web_link(web_links, epid, stacks, sec, poly):
    enormal = sec.GetCellData().GetArray("Normals")
    mat = poly.GetCellData().GetArray("mat")
    angle = poly.GetCellData().GetArray("angle")
    cells = poly.GetPolys()

    for i in web_links:
        # print(i)
        to1 = i[2][0]  # connected element 1
        to2 = i[2][1]  # connected element 2
        s1 = stacks[to1]  # stack of element 1
        s2 = stacks[to2]  # stack of connected element 2
        s0 = stacks[i[0]]  # stack of web element
        n11 = epid[to1] + 2 * len(s1) - 1
        n12 = epid[to1] + 2 * len(s1) - 2
        n21 = epid[to2] + 2 * len(s2) - 1
        n22 = epid[to2] + 2 * len(s2) - 2

        n1, n2 = enormal.GetTuple3(i[0]), enormal.GetTuple3(to2)

        if max(n1[0] * n2[0], n1[1] * n2[1]) > 0.5:
            linkto = [j + epid[to2] for j in range(0, len(s2) * 2, 2)]
        else:
            linkto = [n11, n12, n21, n22]

        nds = linkto + [j + epid[i[0]] + i[3] for j in range(0, len(s0) * 2, 2)]

        pts = vtk.vtkPoints()

        p1 = poly.GetPoint(n11)
        p2 = poly.GetPoint(nds[4])

        if linkto[0] == n11:
            dx = -1 if p1[0] < p2[0] else 1
        else:
            dx = 0

        for j in nds[:4]:
            p = poly.GetPoint(j)
            pts.InsertNextPoint(p[0] + dx, p[1], 0)

        for j in nds[4:]:
            p = poly.GetPoint(j)
            pts.InsertNextPoint(p[0], p[1], 0)

        pol = vtk.vtkPolyData()
        pol.SetPoints(pts)
        d2 = vtk.vtkDelaunay2D()
        d2.SetInputData(pol)

        d2.Update()

        dd = d2.GetOutput()

        for j in range(dd.GetNumberOfCells()):
            c = dd.GetCell(j)
            cells.InsertNextCell(4)
            for k in range(c.GetNumberOfPoints()):
                cells.InsertCellPoint(nds[c.GetPointId(k)])

            cells.InsertCellPoint(nds[c.GetPointId(k)])
            mat.InsertNextTuple1(-1)
            angle.InsertNextTuple1(0)


def create_bondline(y, epid, stacks, sec, poly, bondline_material=20):
    lst = []

    bnds = sec.GetBounds()
    for i in range(sec.GetNumberOfPoints()):
        pnt = sec.GetPoint(i)
        if pnt[1] > bnds[3] - y:
            lst.extend(get_eids(sec, i))

    pts = vtk.vtkPoints()
    nds = []
    for j in set(lst):
        s = stacks[j]
        n1 = epid[j] + 2 * len(s) - 1
        n2 = epid[j] + 2 * len(s) - 2
        pts.InsertNextPoint(poly.GetPoint(n1))
        pts.InsertNextPoint(poly.GetPoint(n2))
        nds.extend((n1, n2))

    pol = vtk.vtkPolyData()
    pol.SetPoints(pts)
    d2 = vtk.vtkDelaunay2D()
    d2.SetInputData(pol)

    dd = d2.GetOutput()

    cells = poly.GetPolys()
    material = poly.GetCellData().GetArray("mat")
    angle = poly.GetCellData().GetArray("angle")

    for j in range(dd.GetNumberOfCells()):
        c = dd.GetCell(j)
        cells.InsertNextCell(4)
        for k in range(c.GetNumberOfPoints()):
            cells.InsertCellPoint(nds[c.GetPointId(k)])
        cells.InsertCellPoint(nds[c.GetPointId(k)])
        material.InsertNextTuple1(bondline_material)
        angle.InsertNextTuple1(0)


def get_local_twist(r, var):
    return np.interp(r, var["z"][1] - np.min(var["z"][1]), var["twist"][1])


def get_local_chord(r, var):
    return np.interp(r, var["z"][1] - np.min(var["z"][1]), var["chord"][1])


def write_vtp(section, vtp):
    wr = vtk.vtkXMLPolyDataWriter()
    wr.SetFileName(vtp)
    wr.SetInputData(section)
    wr.Update()
    wr.Write()
    print(f"writing to {vtp}")


def get_avg_normal(normals):
    mean_normal = np.mean(normals, axis=0)
    mean_normal /= np.linalg.norm(mean_normal)
    return mean_normal


def get_node_avg_normal(connected_cell_normals):
    s = connected_cell_normals.shape
    if s[0] < 3:
        return get_avg_normal(connected_cell_normals), None
    corr = np.zeros((s[0], s[0]))
    for i in range(s[0]):
        for j in range(i + 1, s[0]):
            corr[i, j] = np.abs(
                np.dot(connected_cell_normals[i, :], connected_cell_normals[j, :])
            )
    if corr[0, 1] > np.mean([corr[0, 2], corr[1, 2]]):
        return get_avg_normal(connected_cell_normals[:2, :]), 2
    elif corr[0, 2] > np.mean([corr[0, 1], corr[1, 2]]):
        return get_avg_normal(connected_cell_normals[::2, :]), 1
    else:
        return get_avg_normal(connected_cell_normals[1:, :]), 0


def cut_blade(r, vtu, if_bondline=True, rotz=0, var=None, is2d=False, verbose=False):
    if var is None:
        var = {}
    print("# creating cross section mesh from %s at r=%.3f" % (vtu, r))
    workdir = os.path.dirname(vtu)
    # read in the mesh
    rd = pv.read(vtu)

    local_twist = get_local_twist(r, var)
    local_chord = get_local_chord(r, var)

    sec = rd.slice(normal=[0, 0, 1], origin=[0, 0, r])

    # if not is2d:
    #     # slice the mesh at some point
    #     cutter = vtk.vtkCutter()
    #     plane = vtk.vtkPlane()
    #     plane.SetOrigin(0, 0, r)
    #     plane.SetNormal(0, 0, 1)

    #     cutter.SetCutFunction(plane)
    #     cutter.SetInputData(msh)
    #     cutter.Update()
    #     sec = cutter.GetOutput()
    # else:
    #     sec = msh

    # first we rotate around the Z axis as if rotating the whole blade (compensate for pck angle)
    transform = vtk.vtkTransform()
    transform.RotateZ(rotz)
    transformfilter = vtk.vtkTransformFilter()
    transformfilter.SetTransform(transform)
    transformfilter.SetInputData(sec)
    transformfilter.Update()
    rotated_section = transformfilter.GetOutput()

    pts = np.array(
        [
            rotated_section.GetPoint(i)
            for i in range(rotated_section.GetNumberOfPoints())
        ]
    )
    minxy = np.amin(pts, axis=0)
    maxxy = np.amax(pts, axis=0)
    mid_position = 0.5 * (minxy + maxxy)

    bnds = rotated_section.GetBounds()

    table_out = pd.DataFrame(
        [
            [r]
            + list(mid_position)
            + [local_twist, local_chord]
            + [bnds[1] - bnds[0], bnds[3] - bnds[2]],
        ],
        columns=["r", "xavg", "yavg", "zavg", "twist_angle", "local_chord", "dx", "dy"],
    )

    # get a handle with plydata on points and cells
    ctp = vtk.vtkCellDataToPointData()
    ctp.PassCellDataOn()
    ctp.SetInputData(sec)
    ctp.Update()

    sec = ctp.GetOutput()

    # special case that deals with attachment of added webs that might link
    # crossed to elements (so the slice mesh is not linked up), this code links
    # up the slice by searching in negative y direction until the closest point
    # is no longer the unconnected point
    for i in range(sec.GetNumberOfPoints()):
        pb = sec.GetPoint(i)
        e0 = get_eids(sec, i)
        if len(e0) == 1:
            for j in np.linspace(0, 1, 100):
                p = sec.FindPoint(pb[0], pb[1] - j, pb[2])
                if p != i:
                    sec.GetPoints().SetPoint(i, sec.GetPoint(p))
                    break

    # clean up the mesh
    cln = vtk.vtkCleanPolyData()
    cln.SetInputData(sec)
    cln.Update()

    sec = cln.GetOutput()

    # get the max number of plies

    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    poly = vtk.vtkPolyData()
    poly.SetPoints(points)
    poly.SetPolys(cells)

    material = vtk.vtkIntArray()
    material.SetName("mat")
    poly.GetCellData().AddArray(material)
    angle = vtk.vtkFloatArray()
    angle.SetName("angle")
    poly.GetCellData().AddArray(angle)

    nply = sec.GetCellData().GetArray("n_plies")
    thickness = sec.GetCellData().GetArray("thickness")
    enormal = sec.GetCellData().GetArray("Normals")

    # gather cell-wise ply arrays
    plies = []
    for i in range(sec.GetCellData().GetNumberOfArrays()):
        name = sec.GetCellData().GetArrayName(i)
        if name.startswith("ply"):
            plies.append((name, sec.GetCellData().GetArray(i)))

    pid = 0
    epid = []
    join_nodes = []
    web_links = []
    stacks = []  # stack of plies
    mstacks = []  # material stack
    astacks = []  # angle stack

    # loop over all cells
    for i in range(sec.GetNumberOfCells()):
        # get the id, the nodes attached to the cell, and the elements attached to each node
        epid.append(pid)
        nn = get_nids(sec, i)
        e0 = get_eids(sec, nn[0])
        e1 = get_eids(sec, nn[1])

        e0b, e1b = copy.deepcopy(e0), copy.deepcopy(e1)

        # get the normal in node 0, and if there is a web connected, give back the id in e0 of the web element
        nr0, wn0 = get_node_avg_normal(np.array([enormal.GetTuple3(i) for i in e0]))
        # get the normal in node 1
        nr1, wn1 = get_node_avg_normal(np.array([enormal.GetTuple3(i) for i in e1]))

        if wn0 != None:
            del e0[wn0]
        if wn1 != None:
            del e1[wn1]

        x0, x1 = 0.0, 1.0
        if nply.GetTuple1(e0[0]) != nply.GetTuple1(e0[1]) or thickness.GetTuple1(
            e0[0]
        ) != thickness.GetTuple1(e0[1]):
            join_nodes.append((nn[0], e0))
            x0 = 0.1
        if nply.GetTuple1(e1[0]) != nply.GetTuple1(e1[1]) or thickness.GetTuple1(
            e1[0]
        ) != thickness.GetTuple1(e1[1]):
            x1 = 0.9

        if wn0 != None and i == e0b[wn0]:
            x0, x1 = 0, 1.0
            fact = 1.1
            cell_length = sec.GetCell(i).GetLength2() ** 0.5

            nr0 = np.array(enormal.GetTuple3(i))
            nr0[2] = 0.0
            nr0 /= np.linalg.norm(nr0)
            web_links.append([i, nn[0], e0b, 0])
            for jj in range(3):  # loop over the 3 elements attached to n0
                if wn0 != jj:  # if one isn't the current (web) element
                    th = thickness.GetTuple1(e0b[jj]) * 1e-3  # get the thickness
                    break
            x0 = min(fact * th / cell_length, 0.97)

        if wn1 != None and i == e1b[wn1]:
            x0, x1 = 0, 1.0
            fact = 1.1
            cell_length = sec.GetCell(i).GetLength2() ** 0.5

            nr1 = np.array(enormal.GetTuple3(i))
            nr1[2] = 0.0
            nr1 /= np.linalg.norm(nr1)
            web_links.append([i, nn[1], e1b, 1])
            th = thickness.GetTuple1(e1[0]) * 1e-3
            x1 = max(1.0 - fact * th / cell_length, 0.03)

        p0, p1 = zip(
            *[
                [(1.0 - x0) * ii[0] + x0 * ii[1], (1.0 - x1) * ii[0] + x1 * ii[1]]
                for ii in zip(sec.GetPoint(nn[0]), sec.GetPoint(nn[1]))
            ]
        )

        points.InsertNextPoint(p0), points.InsertNextPoint(p1)
        pid += 2
        stack = [0]
        mstack = []
        astack = []
        for j in plies:
            tup = j[1].GetTuple(i)
            if tup[0] != -1 and tup[1] > 0.01:  # filter out very thin plies (like the
                # ones used for puck postprocessing)
                # if the material and angle are the same as the last, don't put in a new ply,
                # but increase the previous ply thickness, this greatly reduces element count in e.g. sparcaps
                if mstack and (mstack[-1] == int(tup[0]) and astack[-1] == tup[2]):
                    stack[-1] += tup[1]
                    p0 = [k[0] + k[1] * tup[1] * 1e-3 for k in zip(p0, nr0)]
                    p1 = [k[0] + k[1] * tup[1] * 1e-3 for k in zip(p1, nr1)]
                    npts = points.GetNumberOfPoints()
                    points.SetPoint(npts - 2, tuple(p0))
                    points.SetPoint(npts - 1, tuple(p1))
                else:
                    stack.append(tup[1] + stack[-1])
                    mstack.append(int(tup[0]))
                    astack.append(tup[2])
                    p0 = [k[0] + k[1] * tup[1] * 1e-3 for k in zip(p0, nr0)]
                    p1 = [k[0] + k[1] * tup[1] * 1e-3 for k in zip(p1, nr1)]
                    points.InsertNextPoint(tuple(p0))
                    points.InsertNextPoint(tuple(p1))
                    pid += 2

        stacks.append(stack)
        mstacks.append(mstack)
        astacks.append(astack)

    # fill up the "normal" elements
    for i in range(sec.GetNumberOfCells()):
        p1 = epid[i]
        p2 = epid[i] + 1
        for j in range(len(mstacks[i])):
            cells.InsertNextCell(4)
            cells.InsertCellPoint(p1 + 2 * j)
            cells.InsertCellPoint(p1 + 2 * j + 2)
            cells.InsertCellPoint(p2 + 2 * j + 2)
            cells.InsertCellPoint(p2 + 2 * j)
            material.InsertNextTuple1(mstacks[i][j])
            angle.InsertNextTuple1(astacks[i][j])

    # join up the plydrop elements
    join_up(join_nodes, epid, stacks, sec, poly)

    # join up the webs with the shell
    web_link(web_links, epid, stacks, sec, poly)

    # create TE bondline
    if if_bondline:
        create_bondline(0.1, epid, stacks, sec, poly)

    if verbose:
        write_vtp(poly, os.path.join(workdir, "inter_%i.vtp" % (1e3 * r)))

    # clean up the polydata, i.e. remove duplicate points
    cln = vtk.vtkCleanPolyData()
    cln.SetInputData(poly)
    cln.Update()

    #  reapply the transform to the output mesh (run the mesh generator in original orientation)
    transform = vtk.vtkTransform()
    transform.RotateZ(rotz)
    transformfilter = vtk.vtkTransformFilter()
    transformfilter.SetTransform(transform)
    transformfilter.SetInputData(cln.GetOutput())
    transformfilter.Update()
    rotated_section = transformfilter.GetOutput()

    transform = vtk.vtkTransform()
    transform.Translate(-mid_position[0], -mid_position[1], -mid_position[2])
    transformfilter = vtk.vtkTransformFilter()

    transformfilter.SetTransform(transform)
    transformfilter.SetInputData(rotated_section)
    transformfilter.Update()
    zeroed_section = transformfilter.GetOutput()

    transform = vtk.vtkTransform()
    transform.RotateZ(local_twist)
    transformfilter = vtk.vtkTransformFilter()
    transformfilter.SetTransform(transform)

    transformfilter.SetInputData(zeroed_section)
    transformfilter.Update()
    detwisted_section = transformfilter.GetOutput()

    out = detwisted_section
    mat = out.GetCellData().GetArray("mat")
    ang = out.GetCellData().GetArray("angle")

    # get material and angle properties for the triangle elements
    for i in range(out.GetNumberOfCells()):
        c = out.GetCell(i)
        npts = c.GetNumberOfPoints()
        m = mat.GetTuple1(i)
        if npts == 3 and m == -1:
            g = c.GetPointIds()
            nb = []
            for j in [(0, 1), (1, 2), (2, 0)]:
                t = vtk.vtkIdList()
                t.InsertNextId(g.GetId(j[0]))
                t.InsertNextId(g.GetId(j[1]))

                idl = vtk.vtkIdList()
                out.GetCellNeighbors(i, t, idl)

                nb.extend(idl.GetId(k) for k in range(idl.GetNumberOfIds()))
            mat.SetComponent(i, 0, -1)
            ang.SetComponent(i, 0, 0)

            for j in nb:
                cp = out.GetCell(j).GetNumberOfPoints()
                if cp == 4:
                    mat.SetComponent(i, 0, mat.GetTuple1(j))
                    ang.SetComponent(i, 0, ang.GetTuple1(j))

    ang2 = vtk.vtkFloatArray()
    ang2.SetName("angle2")
    ang2.SetNumberOfTuples(out.GetNumberOfCells())
    angle_lookup = {}

    c = 1
    for i in range(out.GetNumberOfCells()):
        cl = out.GetCell(i)

        if cl.GetNumberOfPoints() > 2:
            [pt1, pt2] = [cl.GetEdge(1).GetPoints().GetPoint(j) for j in range(2)]

            vec = np.array([pt2[j] - pt1[j] for j in range(3)])
            vec /= np.linalg.norm(vec)
            dot = vtk.vtkMath.Dot([1, 0, 0], vec)

            angle = math.degrees(math.acos(dot))

            if pt2[1] < pt1[1]:
                angle *= -1

            if cl.GetNumberOfPoints() == 3:  # if
                for jj in range(3):
                    if cl.GetPointId(jj) in angle_lookup:
                        angle = angle_lookup[cl.GetPointId(jj)]
            else:
                angle_lookup[cl.GetPointId(0)] = angle

            ang2.SetComponent(i, 0, angle)

    out.GetCellData().AddArray(ang2)

    if verbose:
        write_vtp(out, os.path.join(workdir, "prerealign_%i.vtp" % (1e3 * r)))

    align_normals(out)

    of = os.path.join(workdir, "msec_%i.vtp" % (1e3 * r))
    wrt = vtk.vtkXMLPolyDataWriter()
    wrt.SetInputData(out)
    wrt.SetFileName(of)
    wrt.Write()
    print(f"# written vtk to {of}")
    table_out.to_csv(
        os.path.join(workdir, "section_location_%i.csv" % (1e3 * r)), index=False
    )

    return of


def cut_blade_parallel(
    vtu, rr, if_bondline, rotz, var, verbose=False, is2d=False, debug=False
):
    var = eval(open(var, "r").read())
    part = partial(
        cut_blade,
        vtu=vtu,
        if_bondline=if_bondline,
        rotz=rotz,
        var=var,
        is2d=is2d,
        verbose=verbose,
    )
    if debug == False:
        pool = multiprocessing.Pool()
        al = pool.map(part, rr)
        pool.close()
        pool.join()
    else:
        al = [part(i) for i in rr]
    return al


def main():
    p = argparse.ArgumentParser(
        description="slice through the 3D FEA mesh (joined) to create 2D solid meshes for computing aeroelastic beam properties"
    )
    p.add_argument("--vtu", help="mesh input")
    p.add_argument("--r", default="[30]", help="radius array")
    p.add_argument("--bondline", default=False)
    p.add_argument("--rotz", type=float, default=0)
    p.add_argument("--var", type=str)
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--debug", action="store_true")
    p.add_argument("--is2d", action="store_true")
    args = p.parse_args()

    rr = eval(args.r)
    cut_blade_parallel(
        args.vtu,
        rr,
        args.bondline,
        args.rotz,
        args.var,
        args.verbose,
        args.is2d,
        args.debug,
    )


if __name__ == "__main__":
    main()
