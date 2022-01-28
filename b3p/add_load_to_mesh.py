#! /usr/bin/env python3

import numpy as np
import vtk
import math
import copy
import pyvista
import yaml
import argparse
import scipy.optimize
from functools import partial
from matplotlib import pyplot


def compute_nodal_forces(nz, target_z, target_moment, fmult=1.0):
    assert fmult ** 2 == 1
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


def main():
    p = argparse.ArgumentParser()
    p.add_argument("yaml")
    p.add_argument("grid")
    args = p.parse_args()

    yml = args.yaml
    gridname = args.grid

    config = yaml.load(open(yml, "r"), Loader=yaml.CLoader)
    grid = pyvista.UnstructuredGrid(gridname)

    fig, ax = pyplot.subplots(2, 1)

    lds = config["loads"]
    for i in lds:
        print("loadcase %s" % i)
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

        ax[0].plot(zmy, bmy, label="my moment backcalc from forces ")
        ax[0].plot(z, my, "o", label="my moment target ")
        ax[0].plot(zmx, bmx, label="mx moment backcalc from forces ")
        ax[0].plot(z, mx, "o", label="mx moment target ")
        ax[0].legend(loc="best")
        ax[1].plot(nz, fx, label="fx, sum=%.2f" % fx.sum())
        ax[1].plot(nz, fy, label="fy, sum=%.2f" % fy.sum())
        ax[1].legend(loc="best")

        force_vector = np.zeros_like(grid.points)
        force_vector[loaded_node_ids, 0] = fx
        force_vector[loaded_node_ids, 1] = fy
        print(i)
        grid.point_data["lc_%s" % i] = force_vector

    print("writing loadcases to grid %s" % gridname)
    grid.save(gridname)
    pyplot.savefig("load_output.png")


if __name__ == "__main__":
    main()