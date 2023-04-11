#! /usr/bin/env python3

import pyvista as pv
import fire
import numpy as np
import json
import os
import yaml
import pandas as pd


def load_mm(matmap):
    """Load the material map and material database from a material map file.

    :param matmap: path to the material map file"""
    mm = json.load(open(matmap))
    mdbpath = os.path.join(os.path.dirname(matmap), mm["matdb"])
    mdb = yaml.load(open(mdbpath, "r"), Loader=yaml.FullLoader)
    print(f"Loaded material database from {mdbpath}")
    return mm, mdb


def drape_summary(vtu, matmap=None):
    """Summarize the drape results from a vtu file

    :param vtu: path to the vtu file
    :param matmap: path to the material map file, if not given, then search for material_map.json in the same directory as the vtu file"""
    gr = pv.read(vtu).compute_cell_sizes()
    pl = [i for i in gr.cell_data.keys() if i.startswith("ply_")]

    mat = np.stack([gr.cell_data[i] for i in pl])

    area = gr.cell_data["Area"]
    radius = gr.cell_data["radius"]

    mm, mdb = load_mm(matmap or os.path.join(os.path.dirname(vtu), "material_map.json"))

    mminv = {v: k for k, v in mm.items()}
    total_volume = 1e-3 * mat[:, :, 1] * area * (mat[:, :, 0] > 0)
    mat_used = np.unique(mat[:, :, 0])

    mass_moment = np.zeros_like(radius)

    out = {}
    for i in mat_used:
        if i > 0:
            material = mminv[int(i)]
            mdbv = mdb[material]
            rho = mdbv["rho"] if "rho" in mdbv else mdbv["density"]
            mat_vol = 1e-3 * mat[:, :, 1] * area * (mat[:, :, 0] == i)
            mat_mass = mat_vol * rho
            mass_moment += mat_mass.sum(axis=0) * radius
            out[material] = {
                "id": int(i),
                "material": material,
                "name": mdbv["name"],
                "rho": rho,
                "volume": mat_vol.sum(),
                "mass": mat_mass.sum(),
            }

    dt = pd.DataFrame(out).T
    print("Drape Summary:")
    print(dt)
    print("Total Volume and Mass:")
    print(dt.sum()[["volume", "mass"]].T)
    print(f"Total volume backcheck: {total_volume.sum():.4f} m^3")
    sum_mass_mom = mass_moment.sum()
    print(f"Mass_moment: {sum_mass_mom} kg*m, radius {sum_mass_mom/dt.sum()['mass']}")
    return dt


if __name__ == "__main__":
    fire.Fire(drape_summary)
