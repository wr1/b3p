#! /usr/bin/env python3

# import vtk
import argparse
import numpy as np
import pyvista as pv
import time
import multiprocessing

# meshes = []


def add_missing_data(inp):
    mesh, pd, cd = inp
    for i in pd:
        mesh.point_data.set_array(i[1], i[0])
    for i in cd:
        mesh.cell_data.set_array(i[1], i[0])
    return mesh


def main():
    # global meshes
    p = argparse.ArgumentParser(
        description="Join a series of meshes, i.e. a shell and n web meshes together into a single vtu"
    )
    p.add_argument("meshes", nargs="*")
    p.add_argument("--out", default="__joined_mesh.vtu", help="output file name")
    args = p.parse_args()

    meshes = []
    for i in args.meshes:
        meshes.append(pv.read(i))

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

    # add the zero arrays in parallel
    pool = multiprocessing.Pool()
    cmeshes = pool.map(add_missing_data, dist)
    toc = time.time()

    out = cmeshes[0].merge(cmeshes[1:])

    toc2 = time.time()

    print("time adding missing arrays: ", toc - tic, "\ntime merging:", toc2 - toc)

    out.save(args.out)
    print("written mesh to %s" % args.out)


if __name__ == "__main__":
    main()
