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

# def cmx(fy, nz, z_load):
#     return np.array([np.sum((nz > i) * (nz - i) * -fy) for i in z_load])


# def cmy(fx, nz, z_load):
#     return np.array([np.sum((nz > i) * (nz - i) * fx) for i in z_load])


# def compute_mx_error(fy, nz, z_load, mx_target, smooth=1e-2):
#     mx = cmx(fy, nz, z_load)
#     mx_error = np.array(mx_target) - mx
#     return np.sum(mx_error ** 2) + smooth * sum(fy ** 2)


# def compute_my_error(fx, nz, z_load, my_target, smooth=1e-2):
#     my = cmy(fx, nz, z_load)
#     my_error = np.array(my_target) - my
#     return np.sum(my_error ** 2) + smooth * sum(fx ** 2)


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

    # yml = "blade_test.yml"
    # gridname = "temp/joinedup_mesh.vtu"

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
        # fx, fy = np.zeros_like(nz), np.zeros_like(nz)
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