#! /usr/bin/env python3

import numpy as np
import pyvista
from matplotlib import pyplot as plt
import pandas as pd


def compute_nodal_forces(nz, target_z, target_moment, fmult=1.0):
    assert fmult**2 == 1
    loaded = np.zeros_like(nz).astype(bool)
    force = np.zeros_like(nz)
    for i in reversed(list(zip(target_z, target_moment))):
        # compute the moment coming from the previously applied forces
        initial_moment = ((nz - i[0]) * fmult * force * loaded).sum()
        # get the nodes that are outboard from this point, but have not had a load applied yet
        rel_nodes = (nz > i[0]) & (loaded == 0)
        # compute distance for these nodes
        dz = (nz - i[0]) * rel_nodes
        # compute force needed to match moment curve
        fx_l = fmult * (i[1] - initial_moment) / dz.sum()
        # apply force
        force[np.where(rel_nodes)] = fx_l
        # register nodes as loaded
        loaded += rel_nodes

    zpl = np.unique(nz)
    moment = [((nz > j) * (nz - j) * force * fmult).sum() for j in zpl]

    return force, zpl, moment


def add_load_to_mesh(config, gridname, plotfile=None):
    grid = pyvista.UnstructuredGrid(gridname)

    if plotfile:
        fig, ax = plt.subplots(2, 1, figsize=(18, 13))

    lds = config["loads"]
    for i in lds:
        print(f"** loadcase {i}")
        # get applicable nodes
        for n, j in enumerate(config["loads"][i]["apply"]):
            key, [mn, mx] = j, config["loads"][i]["apply"][j]
            isloaded_thisrule = (grid.point_data[key] > mn) & (
                grid.point_data[key] < mx
            )
            if n == 0:
                is_loaded_node = isloaded_thisrule
            else:
                is_loaded_node = is_loaded_node & isloaded_thisrule

        loaded_node_ids = np.where(is_loaded_node)

        nz = grid.points[loaded_node_ids][:, 2]
        z = lds[i]["z"]
        mx = lds[i]["mx"]
        my = lds[i]["my"]

        assert len(z) == len(mx) == len(my)

        fx, zmy, bmy = compute_nodal_forces(nz, z, my, fmult=1.0)
        fy, zmx, bmx = compute_nodal_forces(nz, z, mx, fmult=-1.0)

        if plotfile:
            ax[0].plot(zmy, bmy, label=None)
            ax[0].plot(z, my, "o", label=f"my lc {i}")
            ax[0].plot(zmx, bmx, label=None)
            ax[0].plot(z, mx, "o", label=f"mx ls {i}")
            ax[0].legend(loc="best")
            ax[1].plot(nz, fx, label="fx, sum=%.2f" % fx.sum())
            ax[1].plot(nz, fy, label="fy, sum=%.2f" % fy.sum())
            ax[1].legend(loc="best")

        force_vector = np.zeros_like(grid.points)
        force_vector[loaded_node_ids, 0] = fx
        force_vector[loaded_node_ids, 1] = fy
        grid.point_data[f"lc_{i}"] = force_vector

    print(f"writing loadcases to grid {gridname}")
    grid.save(gridname)

    if plotfile:
        fig.savefig(plotfile)
        print(f"** written load plot to {plotfile}")
    return grid


if __name__ == "__main__":
    main()
