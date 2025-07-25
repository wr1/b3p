#! /usr/bin/env python3

import argparse
import pandas as pd
import pickle
import numpy as np
import pyvista
import logging
import multiprocessing
from rich.progress import Progress  # Replace tqdm with rich progress

logger = logging.getLogger(__name__)


def get_slab_cover(inp):
    """Get the cover of the ply"""

    # create a boolean array with n_cell rows and n_ply columns presenting true where the ply covers the cell in the chordwise direction
    name, cover, numbering, rr, stack, df = inp
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
    return names[: data.shape[2]], data


def drape_mesh(vtp, stack, key, output_file):
    """drape a stack of plies onto a mesh"""
    logger.info(f"reading {vtp}")
    x = pyvista.read(vtp)

    o = x.point_data_to_cell_data(pass_point_data=True)
    o = o.compute_cell_sizes(length=False, volume=False)
    o = o.compute_normals(cell_normals=False, flip_normals=True)

    # get all the coordinate systems in a dataframe
    df = pd.DataFrame()
    for i in o.cell_data.keys():
        df[i] = o.cell_data[i]

    logger.info("computing ply coverage")

    # Prepare inputs for parallel processing
    slab_inputs = [
        (i["name"], i["cover"], i["numbering"], i["r"], i["stack"], df)
        for i in stack
        if key.strip() == i["grid"] and i["stack"] != []
    ]

    # Run get_slab_cover in parallel with progress bar
    logger.info("Processing slabs in parallel")
    with multiprocessing.Pool() as pool:
        with Progress() as progress:  # Replace tqdm with rich Progress
            task = progress.add_task("Computing ply coverage", total=len(slab_inputs))
            slab_data = list(
                pool.imap(get_slab_cover, slab_inputs)
            )  # Progress is handled via the context
            progress.update(task, completed=len(slab_data))  # Ensure progress updates

    # Combine results with slab names
    slab_data = [
        [i["name"], names, data]
        for i, (names, data) in zip(
            [s for s in stack if key.strip() == s["grid"] and s["stack"] != []],
            slab_data,
        )
    ]

    logger.info("** assigning ply data to grid")

    # add arrays for each ply, and for each slab
    total_thickness = np.zeros_like(df.radius).astype(np.float32)
    n_plies = np.zeros_like(df.radius).astype(int)
    for i in slab_data:
        slabname, ply_names, dat = i
        for n, j in enumerate(ply_names):
            o.cell_data[j] = dat[:, :, n].T

        s_thick = dat[1, :, :].sum(axis=1)
        o.cell_data[f"slab_thickness_{slabname}"] = s_thick
        total_thickness += s_thick
        n_plies += s_thick > 0.0

    # add an array for every material, showing the thickness of the material
    for i in o.cell_data:
        if i.startswith("ply_"):
            mat, t = o.cell_data[i][:, 0], o.cell_data[i][:, 1]
            mid = np.unique(mat).astype(int)[-1]

            arr = f"mat_{mid}_thickness"
            if arr not in o.cell_data:
                o.cell_data[arr] = np.zeros_like(t)

            o.cell_data[arr] += t

    o.cell_data["n_plies"] = n_plies
    o.cell_data["is_web"] = key.lower() != "shell"
    o.cell_data["thickness"] = total_thickness
    o.compute_normals(cell_normals=True, inplace=True)
    o = o.compute_cell_sizes()
    n = o.cell_data["Normals"]
    z = np.stack([np.zeros(len(n)), np.zeros(len(n)), -np.ones(len(n))]).T
    y = np.cross(n, z)
    x = np.cross(n, y)

    o.cell_data["y_dir"] = y
    o.cell_data["x_dir"] = x

    output_grid = pyvista.UnstructuredGrid(o)
    output_grid.save(output_file, binary=True)
    logger.info(f"written to {output_file}")
    return output_grid


def main():
    p = argparse.ArgumentParser()
    p.add_argument("pck", help="pickle file containing the plies")
    p.add_argument("vtp", help="_web.vtp file")
    p.add_argument("--key", default="blade", help="blade web name on which to drape")
    p.add_argument("--out", default="__draped.vtu", help="output file")
    args = p.parse_args()

    drape_mesh(args.vtp, pickle.load(open(args.pck, "rb")), args.key, args.out)


if __name__ == "__main__":
    main()
