#! /usr/bin/env python

import pyvista as pv
import argparse
from matplotlib import pyplot as plt
import numpy as np


def plot_thresholded(sm, ax):
    p = sm.cell_centers().points
    # get ply arrays
    pk = [i for i in sm.cell_data.keys() if i.startswith("ply_")]
    # get material ids
    mat = [int(sm.cell_data[i][:, 0].max()) for i in pk]
    pt = np.stack(sm.cell_data[i][:, 1] for i in pk)

    # f, ax = plt.subplots(1, 1, figsize=(20, 12))
    base = np.zeros(pt.shape[1])
    for i in range(pt.shape[0]):
        ax.fill_between(
            p[:, 2],
            y1=base,
            y2=base + pt[i, :],
            label=f"{mat[i]}",
            color=f"C{max(mat[i],0)}",
            alpha=0.3,
        )
        base += pt[i, :]

    # get unique labels into legend
    handles, labels = ax.get_legend_handles_labels()
    temp = dict(zip(labels, handles))
    ax.legend(temp.values(), temp.keys(), loc="best")

    # return ax
    # f.savefig("__temp.png")


if __name__ == "__main__":
    # parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", help="path to vtu file")
    parser.add_argument("--out", help="path to output file", default="__temp.png")
    args = parser.parse_args()

    # build the blade mesh
    mesh = pv.read(args.mesh)

    fig, ax = plt.subplots(4, 1, figsize=(40, 30))

    cnt = mesh.contour(
        scalars="d_w1",
        isosurfaces=[0.2],
    )

    plot_thresholded(cnt, ax[0])

    cnt = mesh.contour(
        scalars="d_te",
        isosurfaces=[0.025],
    )
    print(cnt)

    plot_thresholded(cnt, ax[1])

    cnt = mesh.contour(
        scalars="d_w2_r",
        isosurfaces=[0.01],
    )

    plot_thresholded(cnt, ax[2])

    cnt = mesh.contour(
        scalars="d_le_r",
        isosurfaces=[0.001],
    )

    plot_thresholded(cnt, ax[3])

    fig.savefig(args.out)
    print(f"written output to {args.out}")
