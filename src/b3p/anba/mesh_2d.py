#! /usr/bin/env python3

import vtk
import os
import multiprocessing
import argparse
import numpy as np
import pandas as pd
import math
import copy
import json
import pyvista as pv
import logging

logger = logging.getLogger(__name__)


def get_eids(sec, n):
    """Get cell IDs connected to a point."""
    g = vtk.vtkIdList()
    sec.GetPointCells(n, g)
    return [g.GetId(i) for i in range(g.GetNumberOfIds())]


def get_nids(sec, e):
    """Get point IDs for a given cell."""
    l = vtk.vtkIdList()
    sec.GetCellPoints(e, l)
    return [l.GetId(i) for i in range(l.GetNumberOfIds())]


def join_up(nodes, epid, stacks, sec, poly):
    """Link nodes at ply drops using Delaunay triangulation with scaled thicknesses."""
    cells = poly.GetPolys()
    mat = poly.GetCellData().GetArray("mat")
    angle = poly.GetCellData().GetArray("angle")
    for i, el in nodes:
        s0, s1 = (
            np.array(stacks[el[0]], dtype=float),
            np.array(stacks[el[1]], dtype=float),
        )
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

        x0 = 1 if get_nids(sec, el[0])[0] != get_nids(sec, el[1])[1] else 0
        nmap = list(range(epid[el[0]] + x0, epid[el[0]] + x0 + len(s0) * 2, 2)) + list(
            range(epid[el[1]] + 1 - x0, epid[el[1]] + 1 - x0 + len(s1) * 2, 2)
        )
        for j in range(dd.GetNumberOfCells()):
            c = dd.GetCell(j)
            cells.InsertNextCell(4)
            for k in range(c.GetNumberOfPoints()):
                cells.InsertCellPoint(nmap[c.GetPointId(k)])
            cells.InsertCellPoint(nmap[c.GetPointId(k)])
            mat.InsertNextTuple1(-1)
            angle.InsertNextTuple1(0)


def align_normals(align):
    """Reverse cells with normals pointing in the positive z-direction."""
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
    """Link web elements to shell elements."""
    enormal = sec.GetCellData().GetArray("Normals")
    mat = poly.GetCellData().GetArray("mat")
    angle = poly.GetCellData().GetArray("angle")
    cells = poly.GetPolys()

    for i in web_links:
        to1, to2 = i[2][0], i[2][1]
        s1, s2, s0 = stacks[to1], stacks[to2], stacks[i[0]]
        n11, n12 = epid[to1] + 2 * len(s1) - 1, epid[to1] + 2 * len(s1) - 2
        n21, n22 = epid[to2] + 2 * len(s2) - 1, epid[to2] + 2 * len(s2) - 2
        n1, n2 = enormal.GetTuple3(i[0]), enormal.GetTuple3(to2)

        linkto = (
            [j + epid[to2] for j in range(0, len(s2) * 2, 2)]
            if max(n1[0] * n2[0], n1[1] * n2[1]) > 0.5
            else [n11, n12, n21, n22]
        )
        nds = linkto + [j + epid[i[0]] + i[3] for j in range(0, len(s0) * 2, 2)]

        pts = vtk.vtkPoints()
        p1, p2 = poly.GetPoint(n11), poly.GetPoint(nds[4])
        dx = -1 if linkto[0] == n11 and p1[0] < p2[0] else 1 if linkto[0] == n11 else 0

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
    """Create a bondline at the trailing edge."""
    lst = []
    bnds = sec.GetBounds()
    for i in range(sec.GetNumberOfPoints()):
        if sec.GetPoint(i)[1] > bnds[3] - y:
            lst.extend(get_eids(sec, i))

    pts = vtk.vtkPoints()
    nds = []
    for j in set(lst):
        s = stacks[j]
        n1, n2 = epid[j] + 2 * len(s) - 1, epid[j] + 2 * len(s) - 2
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


def get_local_thickness(r, var):
    return np.interp(r, var["z"][1] - np.min(var["z"][1]), var["thickness"][1])


def write_vtp(section, vtp):
    """Write a VTK PolyData object to a .vtp file."""
    wr = vtk.vtkXMLPolyDataWriter()
    wr.SetFileName(vtp)
    wr.SetInputData(section)
    wr.Update()
    wr.Write()
    logger.info(f"writing to {vtp}")


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


def slice_mesh(vtu, r, rotz):
    """Slice the 3D mesh at radius r and apply initial rotation."""
    rd = pv.read(vtu)  # Using PyVista for reading, then converting to VTK
    sec = rd.slice(normal=[0, 0, 1], origin=[0, 0, r])
    transform = vtk.vtkTransform()
    transform.RotateZ(rotz)
    transformfilter = vtk.vtkTransformFilter()
    transformfilter.SetTransform(transform)
    transformfilter.SetInputData(sec)
    transformfilter.Update()
    return transformfilter.GetOutput()


def compute_section_properties(sec, r, var):
    """Compute mid-position and geometric properties of the section."""
    pts = np.array([sec.GetPoint(i) for i in range(sec.GetNumberOfPoints())])
    mid_position = 0.5 * (np.amin(pts, axis=0) + np.amax(pts, axis=0))
    bnds = sec.GetBounds()
    return pd.DataFrame(
        [
            [r]
            + list(mid_position)
            + [
                get_local_twist(r, var),
                get_local_chord(r, var),
                get_local_thickness(r, var),
            ]
            + [bnds[1] - bnds[0], bnds[3] - bnds[2]]
        ],
        columns=[
            "r",
            "xavg",
            "yavg",
            "zavg",
            "twist_angle",
            "local_chord",
            "local_thickness",
            "dx",
            "dy",
        ],
    ), mid_position


def link_unconnected_points(sec):
    """Link unconnected points by adjusting their positions."""
    for i in range(sec.GetNumberOfPoints()):
        pb = sec.GetPoint(i)
        e0 = get_eids(sec, i)
        if len(e0) == 1:
            for j in np.linspace(0, 1, 100):
                p = sec.FindPoint(pb[0], pb[1] - j, pb[2])
                if p != i:
                    sec.GetPoints().SetPoint(i, sec.GetPoint(p))
                    break
    cln = vtk.vtkCleanPolyData()
    cln.SetInputData(sec)
    cln.Update()
    return cln.GetOutput()


def build_2d_mesh(sec):
    """Build the 2D mesh with ply stacks."""
    poly = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
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
    plies = [
        (sec.GetCellData().GetArrayName(i), sec.GetCellData().GetArray(i))
        for i in range(sec.GetCellData().GetNumberOfArrays())
        if sec.GetCellData().GetArrayName(i).startswith("ply")
    ]

    pid = 0
    epid = []
    join_nodes = []
    web_links = []
    stacks, mstacks, astacks = [], [], []

    for i in range(sec.GetNumberOfCells()):
        epid.append(pid)
        nn = get_nids(sec, i)
        e0, e1 = get_eids(sec, nn[0]), get_eids(sec, nn[1])
        e0b, e1b = copy.deepcopy(e0), copy.deepcopy(e1)

        nr0, wn0 = get_node_avg_normal(np.array([enormal.GetTuple3(e) for e in e0]))
        nr1, wn1 = get_node_avg_normal(np.array([enormal.GetTuple3(e) for e in e1]))

        if wn0 is not None:
            del e0[wn0]
        if wn1 is not None:
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

        if wn0 is not None and i == e0b[wn0]:
            x0, x1 = 0, 1.0
            fact = 1.1
            cell_length = sec.GetCell(i).GetLength2() ** 0.5
            nr0 = np.array(enormal.GetTuple3(i))
            nr0[2] = 0.0
            nr0 /= np.linalg.norm(nr0)
            web_links.append([i, nn[0], e0b, 0])
            th = thickness.GetTuple1(e0b[0 if wn0 != 0 else 1]) * 1e-3
            x0 = min(fact * th / cell_length, 0.97)

        if wn1 is not None and i == e1b[wn1]:
            x0, x1 = 0, 1.0
            fact = 1.1
            cell_length = sec.GetCell(i).GetLength2() ** 0.5
            nr1 = np.array(enormal.GetTuple3(i))
            nr1[2] = 0.0
            nr1 /= np.linalg.norm(nr1)
            web_links.append([i, nn[1], e1b, 1])
            th = thickness.GetTuple1(e1[0]) * 1e-3
            x1 = max(1.0 - fact * th / cell_length, 0.03)

        p0 = tuple(
            (1.0 - x0) * a + x0 * b
            for a, b in zip(sec.GetPoint(nn[0]), sec.GetPoint(nn[1]))
        )
        p1 = tuple(
            (1.0 - x1) * a + x1 * b
            for a, b in zip(sec.GetPoint(nn[0]), sec.GetPoint(nn[1]))
        )
        points.InsertNextPoint(p0)
        points.InsertNextPoint(p1)
        pid += 2

        stack, mstack, astack = [0], [], []
        for _, ply_data in plies:
            tup = ply_data.GetTuple(i)
            if tup[0] != -1 and tup[1] > 0.01:
                if mstack and mstack[-1] == int(tup[0]) and astack[-1] == tup[2]:
                    stack[-1] += tup[1]
                    p0 = tuple(k[0] + k[1] * tup[1] * 1e-3 for k in zip(p0, nr0))
                    p1 = tuple(k[0] + k[1] * tup[1] * 1e-3 for k in zip(p1, nr1))
                    npts = points.GetNumberOfPoints()
                    points.SetPoint(npts - 2, p0)
                    points.SetPoint(npts - 1, p1)
                else:
                    stack.append(tup[1] + stack[-1])
                    mstack.append(int(tup[0]))
                    astack.append(tup[2])
                    p0 = tuple(k[0] + k[1] * tup[1] * 1e-3 for k in zip(p0, nr0))
                    p1 = tuple(k[0] + k[1] * tup[1] * 1e-3 for k in zip(p1, nr1))
                    points.InsertNextPoint(p0)
                    points.InsertNextPoint(p1)
                    pid += 2

        stacks.append(stack)
        mstacks.append(mstack)
        astacks.append(astack)

    for i in range(sec.GetNumberOfCells()):
        p1, p2 = epid[i], epid[i] + 1
        for j in range(len(mstacks[i])):
            cells.InsertNextCell(4)
            cells.InsertCellPoint(p1 + 2 * j)
            cells.InsertCellPoint(p1 + 2 * j + 2)
            cells.InsertCellPoint(p2 + 2 * j + 2)
            cells.InsertCellPoint(p2 + 2 * j)
            material.InsertNextTuple1(mstacks[i][j])
            angle.InsertNextTuple1(astacks[i][j])

    return poly, epid, stacks, join_nodes, web_links


# def apply_transforms(poly, rotz, mid_position, local_twist):
#     """Apply transformations to align the 2D mesh in a single step."""
#     cln = vtk.vtkCleanPolyData()
#     cln.SetInputData(poly)
#     cln.Update()

#     transform = vtk.vtkTransform()
#     transform.RotateZ(rotz)
#     transform.Translate(-mid_position[0], -mid_position[1], -mid_position[2])
#     transform.RotateZ(local_twist)
#     transformfilter = vtk.vtkTransformFilter()
#     transformfilter.SetTransform(transform)
#     transformfilter.SetInputData(cln.GetOutput())
#     transformfilter.Update()
#     return transformfilter.GetOutput()


# def apply_transforms(poly, rotz, mid_position, local_twist):
#     """Apply transformations to align the 2D mesh in a single step."""
#     cln = vtk.vtkCleanPolyData()
#     cln.SetInputData(poly)
#     cln.Update()

#     transform = vtk.vtkTransform()
#     transform.PreMultiply()  # Explicitly set to ensure matrix composition order
#     transform.RotateZ(rotz)
#     transform.Translate(-mid_position[0], -mid_position[1], -mid_position[2])
#     transform.RotateZ(local_twist)
#     transformfilter = vtk.vtkTransformFilter()
#     transformfilter.SetTransform(transform)
#     transformfilter.SetInputData(cln.GetOutput())
#     transformfilter.Update()
#     return transformfilter.GetOutput()


def apply_transforms(poly, rotz, mid_position, local_twist):
    """Apply transformations to align the 2D mesh."""
    cln = vtk.vtkCleanPolyData()
    cln.SetInputData(poly)
    cln.Update()

    transform = vtk.vtkTransform()
    transform.RotateZ(rotz)
    transformfilter = vtk.vtkTransformFilter()
    transformfilter.SetTransform(transform)
    transformfilter.SetInputData(cln.GetOutput())
    transformfilter.Update()

    transform = vtk.vtkTransform()
    transform.Translate(-mid_position[0], -mid_position[1], -mid_position[2])
    transformfilter2 = vtk.vtkTransformFilter()
    transformfilter2.SetTransform(transform)
    transformfilter2.SetInputData(transformfilter.GetOutput())
    transformfilter2.Update()

    transform = vtk.vtkTransform()
    transform.RotateZ(local_twist)
    transformfilter3 = vtk.vtkTransformFilter()
    transformfilter3.SetTransform(transform)
    transformfilter3.SetInputData(transformfilter2.GetOutput())
    transformfilter3.Update()
    return transformfilter3.GetOutput()


def assign_triangle_properties(poly):
    """Assign material and angle properties to triangle elements."""
    mat = poly.GetCellData().GetArray("mat")
    ang = poly.GetCellData().GetArray("angle")
    for i in range(poly.GetNumberOfCells()):
        c = poly.GetCell(i)
        if c.GetNumberOfPoints() == 3 and mat.GetTuple1(i) == -1:
            g = c.GetPointIds()
            nb = []
            for j in [(0, 1), (1, 2), (2, 0)]:
                t = vtk.vtkIdList()
                t.InsertNextId(g.GetId(j[0]))
                t.InsertNextId(g.GetId(j[1]))
                idl = vtk.vtkIdList()
                poly.GetCellNeighbors(i, t, idl)
                nb.extend(idl.GetId(k) for k in range(idl.GetNumberOfIds()))
            mat.SetComponent(i, 0, -1)
            ang.SetComponent(i, 0, 0)
            for j in nb:
                if poly.GetCell(j).GetNumberOfPoints() == 4:
                    mat.SetComponent(i, 0, mat.GetTuple1(j))
                    ang.SetComponent(i, 0, ang.GetTuple1(j))


def compute_angle2(poly):
    """Compute and add angle2 array based on cell orientations."""
    ang2 = vtk.vtkFloatArray()
    ang2.SetName("angle2")
    ang2.SetNumberOfTuples(poly.GetNumberOfCells())
    angle_lookup = {}

    for i in range(poly.GetNumberOfCells()):
        cl = poly.GetCell(i)
        if cl.GetNumberOfPoints() > 2:
            pt1, pt2 = [cl.GetEdge(1).GetPoints().GetPoint(j) for j in range(2)]
            vec = np.array([pt2[j] - pt1[j] for j in range(3)])
            vec /= np.linalg.norm(vec)
            dot = vtk.vtkMath.Dot([1, 0, 0], vec)
            angle = math.degrees(math.acos(dot)) * (-1 if pt2[1] < pt1[1] else 1)

            if cl.GetNumberOfPoints() == 3:
                for jj in range(3):
                    if cl.GetPointId(jj) in angle_lookup:
                        angle = angle_lookup[cl.GetPointId(jj)]
                        break
            else:
                angle_lookup[cl.GetPointId(0)] = angle
            ang2.SetComponent(i, 0, angle)

    poly.GetCellData().AddArray(ang2)


def cut_blade(r, vtu, if_bondline=True, rotz=0, var=None, verbose=False):
    """Create a 2D cross-sectional mesh from a 3D VTU file at radius r."""
    var = var or {}
    workdir = os.path.dirname(vtu)
    if not os.path.exists(os.path.join(workdir, "2d")):
        os.makedirs(os.path.join(workdir, "2d"), exist_ok=True)

    output_file = os.path.join(workdir, "2d", f"msec_{int(1e3 * r)}.vtp")

    if os.path.exists(output_file):
        logger.info(f"2d mesh {output_file} already exists - skipping ")
        return output_file

    logger.info(f"# creating cross section mesh from {vtu} at r={r:.3f}")
    sec = slice_mesh(vtu, r, rotz)
    table_out, mid_position = compute_section_properties(sec, r, var)

    ctp = vtk.vtkCellDataToPointData()
    ctp.PassCellDataOn()
    ctp.SetInputData(sec)
    ctp.Update()
    sec = ctp.GetOutput()

    sec = link_unconnected_points(sec)
    poly, epid, stacks, join_nodes, web_links = build_2d_mesh(sec)

    join_up(join_nodes, epid, stacks, sec, poly)
    web_link(web_links, epid, stacks, sec, poly)
    if if_bondline:
        create_bondline(0.1, epid, stacks, sec, poly)

    if verbose:
        write_vtp(poly, os.path.join(workdir, f"inter_{int(1e3 * r)}.vtp"))

    out = apply_transforms(poly, rotz, mid_position, get_local_twist(r, var))
    assign_triangle_properties(out)
    compute_angle2(out)

    if verbose:
        write_vtp(out, os.path.join(workdir, f"prerealign_{int(1e3 * r)}.vtp"))

    align_normals(out)
    write_vtp(out, output_file)
    table_out.to_csv(
        os.path.join(workdir, "2d", f"section_location_{int(1e3 * r)}.csv"), index=False
    )
    logger.info(f"# written vtk to {output_file}")
    return output_file


def run_parallel_tasks(func, items, args_list):
    """Execute a function over a list of items in parallel or sequentially."""
    if len(items) == 0:
        return []
    parallel = args_list[-1]  # Last argument is parallel flag
    if parallel:
        with multiprocessing.Pool() as pool:
            results = pool.starmap(func, [(item, *args_list[:-1]) for item in items])
        return results
    return [func(item, *args_list[:-1]) for item in items]


def cut_blade_parallel(vtu, rr, if_bondline, rotz, var, verbose=False, parallel=False):
    """Process multiple radii in parallel or sequentially."""
    var_data = json.load(open(var, "r"))
    args = [vtu, if_bondline, rotz, var_data, verbose]
    return run_parallel_tasks(cut_blade, rr, args + [parallel])


def main():
    """Command-line interface for slicing 3D meshes into 2D sections."""
    p = argparse.ArgumentParser(
        description="Slice through a 3D FEA mesh (joined) to create 2D solid meshes for computing aeroelastic beam properties"
    )
    p.add_argument("--vtu", required=True, help="Input VTU mesh file")
    p.add_argument("--r", default="[30]", help="Radius array (e.g., '[10, 20, 30]')")
    p.add_argument("--bondline", action="store_true", help="Include bondline")
    p.add_argument("--rotz", type=float, default=0, help="Z-axis rotation angle")
    p.add_argument("--var", type=str, required=True, help="Path to variable file")
    p.add_argument("--verbose", action="store_true", help="Enable verbose output")
    p.add_argument("--parallel", action="store_true", help="Run in parallel")
    args = p.parse_args()

    rr = eval(args.r)
    cut_blade_parallel(
        args.vtu,
        rr,
        args.bondline,
        args.rotz,
        args.var,
        args.verbose,
        args.parallel,
    )


if __name__ == "__main__":
    main()
