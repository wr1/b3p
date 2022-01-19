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
        initial_moment = ((nz - i[0]) * fmult * force * loaded).sum()
        rel_nodes = (nz > i[0]) & (loaded == 0)
        dz = (nz - i[0]) * rel_nodes
        fx_l = (i[1] - initial_moment) / dz.sum()
        force[np.where(rel_nodes)] = fx_l
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

    lds = config["loads"]
    for i in lds:
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

        pyplot.plot(zmy, bmy, label="my moment backcalc from forces ")
        pyplot.plot(z, my, "o", label="my moment target ")

        pyplot.plot(zmx, bmx, label="mx moment backcalc from forces ")
        pyplot.plot(z, mx, "o", label="mx moment target ")

        pyplot.legend(loc="best")

        force_vector = np.zeros_like(grid.points)
        force_vector[loaded_node_ids, 0] = fx
        force_vector[loaded_node_ids, 1] = fy
        grid.point_data["lc_%s" % i] = force_vector

    print("writing loadcases to grid %s" % gridname)
    grid.save(gridname)
    pyplot.savefig("load_output.png")


if __name__ == "__main__":
    main()