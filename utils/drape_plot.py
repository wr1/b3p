#! /usr/bin/env python

import pyvista as pv
from matplotlib import pyplot as plt
import numpy as np
import fire
import os
import json


def plot_thresholded(sm, ax, title, material_map=None):
    p = sm.cell_centers().points
    # get ply arrays
    pk = [i for i in sm.cell_data.keys() if i.startswith("ply_")]
    # get material ids
    mat = [int(sm.cell_data[i][:, 0].max()) for i in pk]

    pt = np.stack([sm.cell_data[i][:, 1] for i in pk])

    base = np.zeros(pt.shape[1])
    for i in range(pt.shape[0]):
        if mat[i] in material_map:
            ax.fill_between(
                p[:, 2],
                y1=base,
                y2=base + pt[i, :],
                label=f"{material_map[mat[i]]}",
                color=f"C{max(mat[i],0)}",
                alpha=0.3,
            )
            base += pt[i, :]

    ax.set_title(title)
    # get unique labels into legend
    handles, labels = ax.get_legend_handles_labels()
    temp = dict(zip(labels, handles))
    ax.legend(temp.values(), temp.keys(), loc="best")


def drape_plot(meshname, output=None):
    """plot the cross section of the laminate plan on (lengthwise) slices of the blade

    :param mesh: path to the mesh in vtu format, needs to have ply_ arrays in cell data
    """
    assert meshname.endswith(".vtu")
    # build the blade mesh
    mesh = pv.read(meshname)

    if output is None:
        output = meshname.replace(".vtu", "_drapeplot.png")

    mm = json.load(open(os.path.join(os.path.dirname(meshname), "material_map.json")))

    imm = {v: k for k, v in mm.items()}
    imm[-1] = "adhesive"
    # if os.path.exists(os.path.join(os.path.dirname(mesh), 'material_map.json')) is False:
    #     os.remove(output)
    # print(mm)

    fig, ax = plt.subplots(4, 1, figsize=(30, 20))

    cnt = mesh.contour(
        scalars="d_w1",
        isosurfaces=[0.0],
    )

    plot_thresholded(cnt, ax[0], title="sparcap", material_map=imm)

    cnt = mesh.contour(
        scalars="d_te",
        isosurfaces=[0.3],
    )

    plot_thresholded(cnt, ax[1], title="trailing edge ud", material_map=imm)

    cnt = mesh.contour(
        scalars="d_w2_r",
        isosurfaces=[0.01],
    )

    plot_thresholded(cnt, ax[2], title="leading edge panel", material_map=imm)

    cnt = mesh.contour(
        scalars="d_w4_r",
        isosurfaces=[-0.05],
    )

    plot_thresholded(cnt, ax[3], title="trailing edge panel", material_map=imm)

    fig.tight_layout()

    fig.savefig(output)
    print(f"written output to {output}")


def main():
    fire.Fire(drape_plot)


if __name__ == "__main__":
    main()
