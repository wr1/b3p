#! /usr/bin/env python3

import pyvista as pv
import fire
import numpy as np
import json
import os
import yaml
import pandas as pd


class drape_summarizer:
    def sum(self, vtu, matmap=None):
        """Summarize the drape results from a vtu file

        :param vtu: path to the vtu file
        :param matmap: path to the material map file, if not given, then search for material_map.json in the same directory as the vtu file"""
        gr = pv.read(vtu).compute_cell_sizes()
        pl = [i for i in gr.cell_data.keys() if i.startswith("ply_")]

        mat = np.stack([gr.cell_data[i] for i in pl])

        area = gr.cell_data["Area"]

        if matmap == None:
            mm = json.load(
                open(os.path.join(os.path.dirname(vtu), "material_map.json"))
            )
            mdb = yaml.load(
                open(os.path.join(os.path.dirname(vtu), mm["matdb"])),
                Loader=yaml.FullLoader,
            )
        else:
            mm = json.load(open(matmap))
            mdb = yaml.load(
                open(os.path.join(os.path.dirname(matmap), mm["matdb"])),
                Loader=yaml.FullLoader,
            )

        mminv = {v: k for k, v in mm.items()}
        total_volume = 1e-3 * mat[:, :, 1] * area * (mat[:, :, 0] > 0)
        mat_used = np.unique(mat[:, :, 0])

        out = {}
        for i in mat_used:
            if i > 0:
                material = mminv[int(i)]
                mdbv = mdb[material]
                # print(material)
                mat_vol = 1e-3 * mat[:, :, 1] * area * (mat[:, :, 0] == i)
                mat_mass = mat_vol * mdbv["rho"] if "rho" in mdbv else mdbv["density"]
                out[material] = {
                    "id": int(i),
                    "material": material,
                    "name": mdbv["name"],
                    "volume": mat_vol.sum(),
                    "mass": mat_mass.sum(),
                }

        dt = pd.DataFrame(out).T
        print("Drape Summary:")
        print(dt)
        print("Total Volume and Mass:")
        print(dt.sum()[["volume", "mass"]].T)
        print(f"Total volume backcheck: {total_volume.sum():.4f} m^3")


if __name__ == "__main__":
    fire.Fire(drape_summarizer())
