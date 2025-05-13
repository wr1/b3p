#! /usr/bin/env python
import pyvista as pv
import pandas as pd
import numpy as np
import os
import glob
import json
import warnings
import logging

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
    """
    Splits a glueline row of solids into 3 rows through thickness adding intermediate points.
    Parameters:
    - fl (vtk.vtkPolyData): The input glueline as a VTK PolyData object.
    Returns:
    - msh (pyvista.UnstructuredGrid): The resulting mesh after splitting the glueline.
    Raises:
    - None
    """

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

    lkp = np.zeros(upr.max() + 1, dtype=int)
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
    file_path, bondline_width=[[0, 0], [0.5, 0.5], [1, 0.1]], bondline_material_id=0
):
    """
    Add bondline to a VTU file.
    Parameters:
    - file_path (str): The path to the VTU file.
    - bondline_width (list, optional): The bondline width values. Default is [[0, 0], [0.5, 0.5], [1, 0.1]].
    Returns:
    None
    """
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
    # added_points = []
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

    msh = split_glueline(msh)

    msh.cell_data["mat"] = bondline_material_id

    add_zero_arrays(msh, mesh)

    out = pv.merge([mesh, msh])
    of = file_path.replace(".vtu", "_bondline.vtu")
    out.save(of)
    logger.info(f"Saved {of}")


def get_bondline_material(d):
    """
    Retrieve bondline material ID and width from the configuration.

    Args:
        d (dict): Blade dictionary containing configuration data.

    Returns:
        tuple: Bondline material ID and bondline width, or (None, None) if not found.
    """
    wd = d["general"]["workdir"]
    mmap = os.path.join(wd, "drape", "material_map.json")

    if os.path.exists(mmap):
        material_map = glob.glob(mmap)

        mm = json.load(open(material_map[0], "r"))

        if "bondline" in d["mesh"]:
            bondline_width = d["mesh"]["bondline"]["width"]

            bondline_material = d["mesh"]["bondline"]["material"]

            bondline_material_id = mm.get(bondline_material, None)

            if bondline_material_id is None:
                warnings.warn(
                    f"Bondline material '{bondline_material}' not found in material map."
                )
                return None, None

            return bondline_material_id, bondline_width
        else:
            warnings.warn(
                "No bondline configuration found in the mesh section. Proceeding without bondline."
            )
            return None, None
    else:
        warnings.warn("Material map file not found. Proceeding without bondline.")
        return None, None


def add_bondline(bladedict, prefix=None):
    vtu = glob.glob(str(prefix) + "*joined.vtu")
    bondline_material_id, bondline_width = get_bondline_material(bladedict)
    if bondline_material_id is not None:
        add_bondline_to_vtu(
            vtu[0],
            bondline_width=bondline_width,
            bondline_material_id=bondline_material_id,
        )
