#! /usr/bin/env python

import pyvista as pv
from matplotlib import pyplot as plt
import numpy as np
import fire


def plot_thresholded(sm, ax):
    p = sm.cell_centers().points
    # get ply arrays
    pk = [i for i in sm.cell_data.keys() if i.startswith("ply_")]
    # get material ids
    mat = [int(sm.cell_data[i][:, 0].max()) for i in pk]
    pt = np.stack(sm.cell_data[i][:, 1] for i in pk)

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


def drape_plot(mesh, output="__temp.png"):
    """plot the cross section of the laminate plan on (lengthwise) slices of the blade

    :param mesh: path to the mesh in vtu format, needs to have ply_ arrays in cell data
    """
    # build the blade mesh
    mesh = pv.read(mesh)

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

    fig.savefig(output)
    print(f"written output to {output}")


def test_drape_plot():
    drape_plot(
        "temp_portable/test_blade_shell.vtu",
        output="data/blade_0001.png",
    )


if __name__ == "__main__":
    fire.Fire(drape_plot)
