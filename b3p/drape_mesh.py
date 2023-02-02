#! /usr/bin/env python3

import argparse
import pandas as pd
import multiprocessing
import pickle
import vtk
import numpy as np
import time
import pyvista


def get_slab_cover(inp):
    # create a boolean array with n_cell rows and n_ply columns presenting true where the ply covers the cell in the chordwise direction
    name, cover, numbering, rr, stack, df = inp
    # d = pd.DataFrame(cover)
    # print(d)
    # exit()
    one = np.ones_like(rr)

    names = ["ply_%.8i_%s" % (i, name) for i in numbering]

    # get the min and max radius for each ply
    material, thickness, rmin, rmax = np.array(stack).T
    rrt = np.array([df.radius.values]).T

    # compute radius coverage
    rcover = (rrt >= rmin) & (rrt < rmax)

    ply_increments = [cover[i][-1] for i in cover]

    if not any(ply_increments):
        # no increment, so the chordwise cover for each ply is the same
        # create one cover vector
        cov = np.ones_like(df.radius, dtype=bool)
        for i in cover:
            start, end, __ = cover[i]
            over_start = np.interp(df.radius, rr, one * start) <= df[i].values
            under_end = np.interp(df.radius, rr, one * end) >= df[i].values
            cov = cov & over_start & under_end

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

    o = x.point_data_to_cell_data(pass_point_data=True)  # centers.GetOutput()

    o = o.compute_cell_sizes(length=False, volume=False)

    o = o.compute_normals(cell_normals=False, flip_normals=True)

    # get all the coordinate systems in a dataframe
    df = pd.DataFrame()
    for i in o.cell_data.keys():
        df[i] = o.cell_data[i]

    stck = pickle.load(open(args.pck, "rb"))

    print("** computing ply coverage")

    slab_data = [
        [i["name"]]
        + list(
            get_slab_cover(
                (
                    i["name"],
                    i["cover"],
                    i["numbering"],
                    i["r"],
                    i["stack"],
                    df,
                )
            )
        )
        for i in stck
        if args.key.strip() == i["grid"] and i["stack"] != []
    ]
    print("** assigning ply data to grid")
    total_thickness = np.zeros_like(df.radius)
    n_plies = np.zeros_like(df.radius).astype(int)
    for i in slab_data:
        slabname, ply_names, dat = i
        # print(slabname, ply_names)
        for n, j in enumerate(ply_names):
            o.cell_data[j] = dat[:, :, n].T

        s_thick = dat[1, :, :].sum(axis=1)
        o.cell_data[f"slab_thickness_{slabname}"] = s_thick
        total_thickness += s_thick
        n_plies += s_thick > 0.0

    o.cell_data["n_plies"] = n_plies
    o.cell_data["is_web"] = args.key.lower().find("web") != -1

    o.cell_data["thickness"] = total_thickness

    o.compute_normals(cell_normals=True, inplace=True)
    o = o.compute_cell_sizes()
    n = o.cell_data["Normals"]
    z = np.stack([np.zeros(len(n)), np.zeros(len(n)), -np.ones(len(n))]).T
    y = np.cross(n, z)
    x = np.cross(n, y)

    o.cell_data["y_dir"] = y
    o.cell_data["x_dir"] = x

    pyvista.UnstructuredGrid(o).save(args.out)
    print(f"** written to {args.out}")


if __name__ == "__main__":
    main()
