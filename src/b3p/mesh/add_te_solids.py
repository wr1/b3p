#! /usr/bin/env python
import pyvista as pv
import pandas as pd
import numpy as np
import glob
import json
import warnings
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def add_zero_arrays(msh, mesh):
    for i in mesh.cell_data.keys():
        s = mesh.cell_data[i].shape
        if i not in msh.cell_data.keys():
            if len(s) == 1:
                msh.cell_data[i] = np.zeros(msh.n_cells, dtype=mesh.cell_data[i].dtype)
            else:
                msh.cell_data[i] = np.zeros(
                    (msh.n_cells, s[1]), dtype=mesh.cell_data[i].dtype
                )
        else:
            msh.cell_data[i] = msh.cell_data[i].astype(mesh.cell_data[i].dtype)

    for i in mesh.point_data.keys():
        s = mesh.point_data[i].shape
        if len(s) == 1:
            msh.point_data[i] = np.zeros(msh.n_points, dtype=mesh.point_data[i].dtype)
        else:
            msh.point_data[i] = np.zeros(
                (msh.n_points, s[1]), dtype=mesh.point_data[i].dtype
            )


def split_glueline(fl):
    """Splits a glueline row of solids into 3 rows through thickness adding intermediate points."""
    logger.debug(fl.n_cells)

    p = fl.points

    cn = fl.cell_connectivity.reshape((fl.n_cells, 8))

    start = cn[:, [1, 0, 5, 4]].flatten()
    end = cn[:, [2, 3, 6, 7]].flatten()

    pr = np.stack([start, end], axis=1)
    upr = np.unique(pr, axis=0)

    ps = p[upr[:, 0]]
    pe = p[upr[:, 1]]

    npp = 2

    spc = np.linspace(0, 1, npp + 2)[1:-1]

    intermediate_points = (
        ps[:, np.newaxis, :] + (pe - ps)[:, np.newaxis, :] * spc[:, np.newaxis]
    )
    added_points = intermediate_points.reshape(-1, 3)

    apids = np.arange(fl.points.shape[0], fl.points.shape[0] + added_points.shape[0])

    lkp = np.zeros(upr.max() + 1, dtype=int)  # .astype(int)  # , dtype=int)
    lkp[upr[:, 0]] = apids[::2]

    end1 = lkp[cn[:, [1, 0, 5, 4]]]

    cn1 = cn.copy()
    cn2 = cn.copy()
    cn[:, [2, 3, 6, 7]] = end1
    cn1[:, [1, 0, 5, 4]] = end1
    cn1[:, [2, 3, 6, 7]] = end1 + 1
    cn2[:, [1, 0, 5, 4]] = end1 + 1

    cells = np.vstack([cn, cn1, cn2])

    n8 = np.array([8 for i in range(cells.shape[0])])

    cells = np.hstack([n8[:, np.newaxis], cells])

    msh = pv.UnstructuredGrid(
        np.array(cells).flatten(),
        np.array([pv.CellType.HEXAHEDRON for i in cells]),
        np.vstack([fl.points, added_points]),
    )
    return msh


def add_bondline_to_vtu(
    file_path: Path,
    bondline_width=[[0, 0], [0.5, 0.5], [1, 0.1]],
    bondline_material_id=0,
) -> Path:
    """Add bondline to a VTU file."""
    logger.info(f"Adding bondline to {file_path} {bondline_width}")
    # Load the VTU file
    mesh = pv.read(file_path)
    mesh.point_data["bondline_width"] = 0.0

    mesh.cell_data["mat"] = -1

    df = pd.DataFrame(mesh.points, columns=["x", "y", "z"])
    df["d_te"] = mesh.point_data["d_te"]
    df["is_web"] = mesh.point_data["is_web"]
    df["d_abs_dist_from_te"] = mesh.point_data["d_abs_dist_from_te"]
    df["rr"] = mesh.point_data["rr"]

    bw = np.array(bondline_width)

    df["bw"] = np.interp(df.rr, bw[:, 0], bw[:, 1])

    shell_pts = df[df.is_web == 0]
    shell_pts = shell_pts.sort_values("d_abs_dist_from_te")
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
        np.array(cells).astype(int).flatten(),
        np.array([pv.CellType.HEXAHEDRON for i in cells]),
        mesh.points,
    )

    msh = split_glueline(msh)

    msh.cell_data["mat"] = bondline_material_id

    add_zero_arrays(msh, mesh)

    out = pv.merge([mesh, msh])
    of = file_path.replace(".vtu", "_bondline.vtu")
    out.save(of)
    logger.info(f"Saved {of}")
    return Path(of)


def get_bondline_material(material_map: Path, bondline_config: dict) -> tuple:
    """
    Retrieve bondline material ID and width from the configuration.
    """
    if not material_map.exists():
        warnings.warn(
            f"Material map file {material_map} does not exist. Proceeding without bondline."
        )
        return None, None

    mm1 = json.load(material_map.open("r"))
    mm = mm1["map"]

    bondline_width = bondline_config["width"]

    bondline_material = bondline_config["material"]

    bondline_material_id = mm.get(bondline_material, None)

    if bondline_material_id is None:
        warnings.warn(
            f"Bondline material '{bondline_material}' not found in material map."
        )
        return None, None

    logger.info(
        f"Bondline material ID: {bondline_material_id}, Width: {bondline_width}"
    )

    return bondline_material_id, bondline_width


def add_bondline(
    vtu_mesh: Path, material_map: Path, bondline_config: dict
):  # bladedict, prefix=None):
    # vtu = glob.glob(str(prefix) + "*joined.vtu")

    bondline_material_id, bondline_width = get_bondline_material(
        material_map, bondline_config=bondline_config
    )
    if bondline_material_id is not None:
        add_bondline_to_vtu(
            vtu_mesh,
            bondline_width=bondline_width,
            bondline_material_id=bondline_material_id,
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Add bondline to VTU file.")

    parser.add_argument("vtu_file", type=str, help="Path to the VTU file.")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    add_bondline_to_vtu(
        args.vtu_file,
    )
