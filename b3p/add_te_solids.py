#! /usr/bin/env python
import fire
import pyvista as pv
import pandas as pd
import numpy as np
from ruamel import yaml
import os
import glob


def add_zero_arrays(msh, mesh):
    for i in mesh.cell_data.keys():
        s = mesh.cell_data[i].shape
        if len(s) == 1:
            msh.cell_data[i] = np.zeros(msh.n_cells, dtype=mesh.cell_data[i].dtype)
        else:
            msh.cell_data[i] = np.zeros(
                (msh.n_cells, s[1]), dtype=mesh.cell_data[i].dtype
            )
    for i in mesh.point_data.keys():
        s = mesh.point_data[i].shape
        if len(s) == 1:
            msh.point_data[i] = np.zeros(msh.n_points, dtype=mesh.point_data[i].dtype)
        else:
            msh.point_data[i] = np.zeros(
                (msh.n_points, s[1]), dtype=mesh.point_data[i].dtype
            )


def add_bondline_to_vtu(file_path, bondline_width=[[0, 0], [0.5, 0.5], [1, 0.1]]):
    # Load the VTU file
    mesh = pv.read(file_path)
    mesh.point_data["bondline_width"] = 0.0
    # bondline_width = 0.4

    df = pd.DataFrame(mesh.points, columns=["x", "y", "z"])
    df["d_te"] = mesh.point_data["d_te"]
    df["is_web"] = mesh.point_data["is_web"]
    df["d_abs_dist_from_te"] = mesh.point_data["d_abs_dist_from_te"]
    df["rr"] = mesh.point_data["rr"]
    # / mesh.point_data["radius"].max()

    bw = np.array(bondline_width)

    df["bw"] = np.interp(df.rr, bw[:, 0], bw[:, 1])

    # print(df.bw)

    shell_pts = df[df.is_web == 0]
    shell_pts.sort_values("d_abs_dist_from_te", inplace=True)
    grz = [g for g in shell_pts.groupby("z")]

    cells = []
    for cg, ng in zip(grz, grz[1:]):
        cgi, ngi = cg[1].index, ng[1].index
        for i in range(200):
            if cg[1].iloc[1 + i]["d_abs_dist_from_te"] > cg[1].bw.max():
                mesh.point_data["bondline_width"][
                    cg[1].index.tolist() + ng[1].index.tolist()
                ] = cg[1].iloc[i]["d_abs_dist_from_te"]

                break
            cells.append(
                [
                    8,
                    cgi[-1 - i],
                    cgi[-2 - i],
                    cgi[1 + i],
                    cgi[0 + i],
                    ngi[-1 - i],
                    ngi[-2 - i],
                    ngi[1 + i],
                    ngi[0 + i],
                ]
            )
    msh = pv.UnstructuredGrid(
        np.array(cells).flatten(),
        np.array([pv.CellType.HEXAHEDRON for i in cells]),
        mesh.points,
    )
    add_zero_arrays(msh, mesh)
    out = pv.merge([mesh, msh])
    of = file_path.replace(".vtu", "_bondline.vtu")
    out.save(of)
    print(f"Saved {of}")


def add_bondline(yml):
    y = yaml.YAML()
    d = y.load(open(yml, "r"))
    vtu = glob.glob(os.path.join(d["general"]["workdir"] + "_portable", "*joined.vtu"))
    bondline_width = d["mesh"]["bondline_width"]
    add_bondline_to_vtu(vtu[0], bondline_width=bondline_width)


def main():
    fire.Fire(add_bondline)


if __name__ == "__main__":
    main()
