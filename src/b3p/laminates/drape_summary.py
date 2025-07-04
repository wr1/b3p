#! /usr/bin/env python3

import pyvista as pv
import numpy as np
import json
import os

# import yaml
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def load_mm(matmap):
    """Load the material map and material database from a material map file."""
    mm = json.load(open(matmap))

    logger.info(f"Loaded material database from {matmap}")
    return mm["map"], mm["matdb"]


def drape_summary(vtu, matmap=None):
    """Summarize the drape results from a vtu file, including ply and slab-based summaries"""
    gr = pv.read(vtu).compute_cell_sizes()
    pl = [i for i in gr.cell_data.keys() if i.startswith("ply_")]
    sl = [i for i in gr.cell_data.keys() if i.startswith("slab_thickness_")]

    mat = np.stack([gr.cell_data[i] for i in pl])
    np.stack([gr.cell_data[i] for i in sl])

    area = gr.cell_data["Area"]
    radius = gr.cell_data["radius"]

    mm, mdb = load_mm(matmap or os.path.join(os.path.dirname(vtu), "material_map.json"))

    mminv = {v: k for k, v in mm.items()}

    # Ply-based summary
    logger.info("Computing ply-based summary")
    ply_volumes = (1e-3 * mat[:, :, 1] * area).sum(axis=1)

    matid = mat[:, :, 0].max(axis=1).astype(int).tolist()
    matmdb = [(mminv[i] if i in mminv else "no mat") for i in matid]
    rho = np.array([(mdb[i]["rho"] if i in mdb else 0) for i in matmdb])

    df_ply = pd.DataFrame(
        zip(pl, ply_volumes, matid, matmdb, rho * ply_volumes),
        columns=["ply", "volume", "matfea", "matmdb", "ply_mass"],
    )
    dirname = os.path.dirname(vtu)
    df_ply.to_csv(os.path.join(dirname, "ply_bom.csv"))

    total_mass_ply = df_ply.ply_mass.sum()

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
    logger.info("Total Volume and Mass (Ply-based):")
    logger.info(dt.sum()[["volume", "mass"]].T)
    logger.info(f"Total volume backcheck (ply-based): {total_volume.sum():.4f} m^3")
    logger.info(f"Total mass backcheck (ply-based): {total_mass_ply:.4f} kg")
    sum_mass_mom = mass_moment.sum()
    logger.info(
        f"Mass_moment: {sum_mass_mom} kg*m, radius {sum_mass_mom / dt.sum()['mass']}"
    )
    return dt
