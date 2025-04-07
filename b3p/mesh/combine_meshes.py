#! /usr/bin/env python3
import argparse
import numpy as np
import pyvista as pv
import time


def add_missing_data(inp):
    mesh, pd, cd = inp
    for i in pd:
        mesh.point_data.set_array(i[1], i[0])
    for i in cd:
        mesh.cell_data.set_array(i[1], i[0])
    return mesh


def is_nonzero_array(arr):
    # check if there is any nonzero entry in the 1 column (thickness)
    if len(arr.shape) == 2 and arr.shape[1] == 3:
        return np.count_nonzero(arr[:, 1]) > 0
    return True


def combine_meshes(meshes, output_filename):
    """Combine a series of meshes into a single vtu file.

    meshes: list of vtu files
    output_filename: name of output file
    """
    meshes = [pv.read(i) for i in meshes]
    # find all the point and cell data arrays
    all_pd = [
        (j, x.point_data[j].shape, x.point_data[j].dtype)
        for x in meshes
        for j in x.point_data.keys()
    ]
    all_cd = [
        (j, x.cell_data[j].shape, x.cell_data[j].dtype)
        for x in meshes
        for j in x.cell_data.keys()
    ]

    tic = time.time()

    # for each mesh, find the missing point and cell arrays and create zero arrays
    dist = []
    for m in meshes:
        da = [m, [], []]
        for j in all_pd:
            if j[0] not in m.point_data:
                a = np.zeros((m.n_points, j[1][1] if len(j[1]) > 1 else 1), dtype=j[2])
                da[1].append((j[0], a))
        for j in all_cd:
            if j[0] not in m.cell_data:
                a = np.zeros((m.n_cells, j[1][1] if len(j[1]) > 1 else 1), dtype=j[2])
                da[2].append((j[0], a))
        dist.append(da)

    toc0 = time.time()

    # this is faster than using multiprocessing,
    # prob because serializing the meshes is slow
    cmeshes = [add_missing_data(i) for i in dist]

    toc1 = time.time()

    out = (
        cmeshes[0]
        .merge(cmeshes[1:])
        .compute_cell_quality(quality_measure="aspect_ratio")
    )

    toc2 = time.time()

    print(
        f"timing: \n\tcreate 0 arrays: { toc0 - tic:.3f} \n\ttime adding: {toc1-toc0:.3f}\n\ttime merging: {toc2-toc1:.3f}"
    )

    out.save(output_filename)
    print(f"written mesh to {output_filename}")
    return out


def main():
    """Combine a series of meshes into a single vtu file."""
    # global meshes
    p = argparse.ArgumentParser(
        description="Join a series of meshes, i.e. a shell and n web meshes together into a single vtu"
    )
    p.add_argument("meshes", nargs="*", help="Input meshes.")
    p.add_argument("--out", default="__joined_mesh.vtu", help="output file name")
    args = p.parse_args()

    combine_meshes(args.meshes, args.out)


if __name__ == "__main__":
    main()
