#! /usr/bin/env python3

import argparse

# import vedo
import pandas as pd
import multiprocessing
import pickle
import vtk
import numpy as np
import time
import pyvista


# def get_ply_cover(inp):
#     name, cover, ply, df = inp
#     rstart, rend = ply[1][1][2:4]
#     ply_key = ply[1][0]
#     material, thickness = ply[1][1][:2]

#     sel = (df.radius >= rstart) & (df.radius <= rend)
#     for k in cover:
#         # for each cover rule, make a boolean array on whether the ply is between max and min for the current rule
#         coord, start, end, inc, r = k[0], k[1], k[2], k[3], k[4]
#         cov = (np.interp(df.radius, r, start + ply[0] * inc) <= df[coord].values) & (
#             np.interp(df.radius, r, end + ply[0] * inc) >= df[coord].values
#         )
#         # the total coverage comes from the radius coverage & each cover rule
#         sel = sel & cov

#     mat = material * sel
#     out = np.stack(
#         [np.where(mat == 0, -1, mat), thickness * sel, np.zeros(len(sel))]
#     ).T.astype(np.float32)

#     return ("ply_%.8i_%s" % (ply_key, name), material, out)


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


def get_slab_cover(inp):
    # create a boolean array with n_cell rows and n_ply columns presenting true where the ply covers the cell in the chordwise direction
    name, cover, numbering, rr, stack, df = inp
    d = pd.DataFrame(cover)
    one = np.ones_like(rr)

    names = ["ply_%.8i_%s" % (i, name) for i in numbering]

    # get the min and max radius for each ply
    material, thickness, rmin, rmax = np.array(stack).T
    rrt = np.array([df.radius.values]).T

    # compute radius coverage
    rcover = (rrt >= rmin) & (rrt <= rmax)

    if d.abs().max(axis=1)[2] == 0:
        # no increment, so the chordwise cover for each ply is the same
        # create one cover vector
        cov = np.ones_like(df.radius, dtype=bool)
        for i in cover:
            start, end, __ = cover[i]
            cov = (
                cov
                & (np.interp(df.radius, rr, one * start) <= df[i].values)
                & (np.interp(df.radius, rr, one * end) >= df[i].values)
            )

        rccov = np.array([cov]).T & rcover
    else:
        # there is an offset for each ply, so create a cover array for each ply
        cov = np.full((len(df.radius), len(numbering)), True)
        for j, __ in enumerate(numbering):
            for i in cover:
                start, end, inc = cover[i]
                cov[:, j] = (
                    cov[:, j]
                    & (np.interp(df.radius, rr, one * start + j * inc) <= df[i].values)
                    & (np.interp(df.radius, rr, one * end + j * inc) >= df[i].values)
                )
        rccov = cov & rcover

    data = np.stack(
        [
            np.where(rccov, material, -1),
            np.where(rccov, thickness, 0),
            np.zeros_like(rccov),
        ]
    ).astype(np.float32)

    return names, data


def main():
    p = argparse.ArgumentParser()
    p.add_argument("pck")
    p.add_argument("vtp", help="_web.vtp file")
    p.add_argument("--key", default="blade")
    p.add_argument("--out", default="__draped.vtu")
    args = p.parse_args()

    x = pyvista.read(args.vtp)

    o = x.point_data_to_cell_data()  # centers.GetOutput()

    o = o.compute_cell_sizes(length=False, volume=False)

    o = o.compute_normals(cell_normals=False, flip_normals=True)

    # get all the coordinate systems in a dataframe
    df = pd.DataFrame()
    for i in o.cell_data.keys():
        df[i] = o.cell_data[i]

    stck = pickle.load(open(args.pck, "rb"))

    print("** computing ply coverage")

    slab_data = []
    for i in stck:  # for each slab
        if args.key.strip() == i["grid"]:  # if the key corresponds to this grid key
            slab_data.append(
                [i["name"]]
                + list(
                    get_slab_cover(
                        (i["name"], i["cover"], i["numbering"], i["r"], i["stack"], df)
                    )
                )
            )

    print("** assigning ply data to grid")
    total_thickness = np.zeros_like(df.radius)
    for i in slab_data:
        slabname, ply_names, dat = i
        for n, j in enumerate(ply_names):
            o.cell_data[j] = dat[:, :, n].T

        s_thick = dat[1, :, :].sum(axis=1)
        o.cell_data["slab_thickness_%s" % slabname] = s_thick
        total_thickness += s_thick

    o.cell_data["is_web"] = args.key.lower().find("web") != -1

    o.cell_data["thickness"] = total_thickness

    print("** writing to file")

    o.save("gaai.vtp")

    exit()

    lst = []
    for i in stck:  # for each slab stack
        # name, grid, cover, stack_numbering, stack = i
        if args.key.strip() == i["grid"]:  # if the key corresponds to this grid's key

            for j in enumerate(
                zip(stack_numbering, stack)
            ):  # for all plies in the slab
                lst.append((name, cover, j, df))

    # only for the incremented plies (i.e. plies that are offset between each ply in the stack)
    # do you need to recalculate cover, so it's faster to precalculate cover

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
        o.cell_data[name] = data
        if material not in mtarr:  # create a thickness array for each material
            mtarr[material] = data[:, 1].astype(np.float32)
        else:
            mtarr[material] += data[:, 1].astype(np.float32)
        n_plies += data[:, 1] > 0

    o.cell_data["n_plies"] = n_plies

    thickness = np.zeros(
        len(data[:, 1]), dtype=np.float32
    )  # add a total thickness array
    for i in mtarr:
        o.cell_data["mat_%i_thickness" % i] = mtarr[i]
        thickness += mtarr[i]

    o.cell_data["thickness"] = thickness  # addCellArray(thickness, "thickness")

    """
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

    o.addCellArray(y, "y_dir")
    o.addCellArray(x, "x_dir")

    # add cell area array
    o.addCellArray(get_area_array(o), "area")

    if args.key.lower().find("web") == -1:
        o.addCellArray(np.zeros(len(n)), "is_web")
    else:
        o.addCellArray(np.ones(len(n)), "is_web")

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
    """


if __name__ == "__main__":
    main()
