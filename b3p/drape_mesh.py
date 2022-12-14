#! /usr/bin/env python3

import argparse
import vedo
import pandas as pd
import multiprocessing
import pickle
import vtk
import numpy as np
import time


def get_ply_cover(inp):
    name, cover, ply, df = inp
    rstart, rend = ply[1][1][2:4]
    ply_key = ply[1][0]
    material, thickness = ply[1][1][:2]

    sel = (df.radius >= rstart) & (df.radius <= rend)
    for k in cover:
        # for each cover rule, make a boolean array on whether the ply is between max and min for the current rule
        coord, start, end, inc, r = k[0], k[1], k[2], k[3], k[4]
        cov = (np.interp(df.radius, r, start + ply[0] * inc) <= df[coord].values) & (
            np.interp(df.radius, r, end + ply[0] * inc) >= df[coord].values
        )
        # the total coverage comes from the radius coverage & each cover rule
        sel = sel & cov

    mat = material * sel
    out = np.stack(
        [np.where(mat == 0, -1, mat), thickness * sel, np.zeros(len(sel))]
    ).T.astype(np.float32)

    return ("ply_%.8i_%s" % (ply_key, name), material, out)


def get_area_array(grid):
    "get array for cell area"
    m = grid.polydata()
    cell_area = []
    for i in range(m.GetNumberOfCells()):
        ids = m.GetCell(i).GetPointIds()
        a1 = vtk.vtkTriangle.TriangleArea(
            m.GetPoint(ids.GetId(0)), m.GetPoint(ids.GetId(1)), m.GetPoint(ids.GetId(2))
        )
        a2 = vtk.vtkTriangle.TriangleArea(
            m.GetPoint(ids.GetId(0)), m.GetPoint(ids.GetId(2)), m.GetPoint(ids.GetId(3))
        )
        cell_area.append(a1 + a2)
    return np.array(cell_area)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("pck")
    p.add_argument("vtp", help="_web.vtp file")
    p.add_argument("--key", default="blade")
    p.add_argument("--out", default="__draped.vtu")
    args = p.parse_args()

    x = vedo.load(args.vtp)

    # translate the coordinates from nodes to cell centers
    centers = vtk.vtkPointDataToCellData()
    centers.SetInputData(x.polydata())
    centers.PassPointDataOn()
    centers.Update()
    poly = centers.GetOutput()

    o = vedo.Mesh(poly)

    df = pd.DataFrame()
    for i in o.celldata.keys():
        df[i] = o.celldata[i]

    stck = pickle.load(open(args.pck, "rb"))

    lst = []
    for i in stck:  # make a list of datasets for each ply
        name, grid, cover, stack_numbering, stack = i
        if args.key.strip() == grid:
            for j in enumerate(zip(stack_numbering, stack)):
                lst.append((name, cover, j, df))

    p = multiprocessing.Pool()  # create the coverage arrays in parallel
    dsets = p.map(get_ply_cover, lst)
    p.close()

    print("# compute ply cell coverage")

    # add the coverage arrays to the mesh
    print("# add ply arrays to grid")
    mtarr = {}
    n_plies = np.zeros(len(dsets[0][2][:, 0]), dtype=int)
    for i in dsets:
        name, material, data = i
        o.celldata[name] = data
        if material not in mtarr:  # create a thickness array for each material
            mtarr[material] = data[:, 1].astype(np.float32)
        else:
            mtarr[material] += data[:, 1].astype(np.float32)
        n_plies += data[:, 1] > 0

    o.celldata["n_plies"] = n_plies

    thickness = np.zeros(
        len(data[:, 1]), dtype=np.float32
    )  # add a total thickness array
    for i in mtarr:
        o.celldata["mat_%i_thickness" % i] = mtarr[i]
        thickness += mtarr[i]

    o.celldata["thickness"] = thickness  # addCellArray(thickness, "thickness")

    o.reverse()

    pdnorm = vtk.vtkPolyDataNormals()
    pdnorm.SetInputData(o.polydata())
    pdnorm.ComputePointNormalsOn()
    pdnorm.ComputeCellNormalsOn()
    pdnorm.FlipNormalsOff()
    pdnorm.ConsistencyOn()
    pdnorm.Update()
    o._update(pdnorm.GetOutput())

    # compute x and y fiber directions
    n = o.normals(cells=True)
    z = np.stack([np.zeros(len(n)), np.zeros(len(n)), -np.ones(len(n))]).T
    y = np.cross(n, z)
    x = np.cross(n, y)

    o.celldata["y_dir"] = y  # addCellArray(y, "y_dir")
    o.celldata["x_dir"] = x  # addCellArray(x, "x_dir")

    # add cell area array
    o.celldata["area"] = get_area_array(o)  # addCellArray(get_area_array(o), "area")

    if args.key.lower().find("web") == -1:
        o.celldata["is_web"] = np.zeros(
            len(n)
        )  # addCellArray(np.zeros(len(n)), "is_web")
    else:
        o.celldata["is_web"] = np.ones(
            len(n)
        )  # addCellArray(np.ones(len(n)), "is_web")

    # translate to unstructuredgrid
    tous = vtk.vtkAppendFilter()
    tous.SetInputData(o.polydata())
    tous.Update()

    # write to vtu
    wr = vtk.vtkXMLUnstructuredGridWriter()
    wr.SetFileName(args.out)
    wr.SetInputData(tous.GetOutput())
    wr.Update()

    print("written to %s" % args.out)


if __name__ == "__main__":
    main()
