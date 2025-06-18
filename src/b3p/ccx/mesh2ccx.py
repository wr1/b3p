#! /usr/bin/env python3

import pyvista as pv
import numpy as np
import vtk
import time
import json
import os

# import yaml
import pandas as pd
import logging
# from rich.logging import RichHandler

logger = logging.getLogger(__name__)
# logger.handlers.clear()
# logger.addHandler(RichHandler(rich_tracebacks=True))


def zero_midside_loads(mesh):
    if mesh.celltypes[0] == 23:
        conn = mesh.cell_connectivity.reshape(
            (
                mesh.GetNumberOfCells(),
                int(mesh.cell_connectivity.shape[0] / mesh.GetNumberOfCells()),
            )
        )
        midsides = conn[:, 4:8].flatten()
        for i in mesh.point_data:
            if i.startswith("lc_"):
                mesh.point_data[i][midsides] *= 0.0
    return mesh


def make_shell_section(elem_id, plyarray, merge_adjacent_plies=True, zero_angle=True):
    plies = []
    filtered_plyarray = plyarray[plyarray[:, 1] > 1e-6]

    for j in filtered_plyarray:
        if plies and plies[-1][1] == j[0] and merge_adjacent_plies:
            plies[-1][0] += j[1] * 1e-3
        else:
            plies.append([j[1] * 1e-3, j[0]])

    if zero_angle:
        section_string = "".join("%f,,m%i,0\n" % tuple(i) for i in plies)
    else:
        section_string = "".join(
            "%f,,m%i,or%i\n" % tuple(i + [elem_id + 1]) for i in plies
        )

    return len(plies), section_string


def material_db_to_ccx(materials, matmap=None, force_iso=False):
    """find the material db and write properties to a ccx block"""
    mat_db = None
    if os.path.isfile(matmap):  # check if the material map file is there
        mm1 = json.load(open(matmap, "r"))
        mm = mm1["map"]
        mat_db = mm1["matdb"]
    else:
        exit("no material map defined")

    mm_inv = {v: k for k, v in mm.items()}

    matblock = ""
    for i in materials:
        if i > 1e-6:
            material_properties = mat_db[mm_inv[int(i)]]
            matblock += (
                f"** material: {mm_inv[int(i)]} {i} {material_properties['name']}\n"
            )

            if "C" in material_properties and not force_iso:
                logger.info(
                    f"{material_properties['name']} is assumed to be orthotropic"
                )
                C = np.array(material_properties["C"])
                matblock += "** orthotropic material\n"
                matblock += "*material,name=m%i\n*elastic,type=ortho\n" % i
                D = C
                D[0, 3] = C[0, 5]
                D[0, 5] = C[0, 3]
                D[1, 3] = C[1, 5]
                D[1, 5] = C[1, 3]
                D[2, 3] = C[2, 5]
                D[2, 5] = C[2, 3]
                D[3, 3] = C[5, 5]
                D[5, 5] = C[3, 3]
                matblock += (
                    f"{D[0, 0]:.4g},{D[0, 1]:.4g},{D[1, 1]:.4g},"
                    + f"{D[0, 2]:.4g},{D[1, 2]:.4g},{D[2, 2]:.4g},"
                    + f"{D[3, 3]:.4g},{D[4, 4]:.4g},\n"
                    + f"{D[5, 5]:.4g},293\n"
                )
            elif "e11" in material_properties and not force_iso:
                logger.info(f"{material_properties['name']} has engineering constants")
                matblock += "** orthotropic material\n"
                matblock += (
                    "*material,name=m%i\n*elastic,type=engineering constants\n" % i
                )
                matblock += (
                    f"{material_properties['e11']:.4g},{material_properties['e22']:.4g},{material_properties['e33']:.4g},"
                    + f"{material_properties['nu12']:.4g},{material_properties['nu13']:.4g},{material_properties['nu23']:.4g},"
                    + f"{material_properties['g12']:.4g},{material_properties['g13']:.4g},\n"
                    + f"{material_properties['g23']:.4g},293\n"
                )
            else:
                logger.info(f"{material_properties['name']} is assumed to be isotropic")
                nu = min(
                    0.45,
                    max(
                        0.1,
                        (
                            float(material_properties["nu"])
                            if "nu" in material_properties
                            else material_properties["nu12"]
                        ),
                    ),
                )
                E = float(
                    material_properties["Ex"]
                    if "Ex" in material_properties
                    else material_properties["E"]
                )
                matblock += "** isotropic material\n"
                matblock += "*material,name=m%i\n*elastic,type=iso\n" % i
                matblock += f"{E:.4g},{nu:.4g},293\n"

    return matblock


def format_eset(name, eids):
    out = f"*elset,elset={name}\n"
    for i in range(len(eids)):
        out += f"{eids[i]}"
        out += "\n" if (i % 16 == 15) else ","
    if out[-1] == ",":
        out = out[:-1] + "\n"
    return out


def compute_ply_groups(grid, prefix):
    gr = ""
    out = []
    n = 1
    for i in grid.cell_data:
        if i.startswith(prefix):
            thickness = grid.cell_data[i][:, 1].max()
            material = grid.cell_data[i][:, 0].max()
            out.append(
                {
                    "Ply Name": i,
                    "Ply ID": n,
                    "Card Image": "PLY",
                    "Mat Name": f"m{int(material)}",
                    "Thickness": thickness,
                    "Orientation": 0,
                    "Output Results": "yes",
                    "TMANU": "",
                    "DRAPE_ID": 0,
                    "ESID": n,
                }
            )
            eids = np.where(grid.cell_data[i][:, 1] > 0)[0] + 1
            gr += format_eset(i, eids)
            n += 1
    return gr, pd.DataFrame(out)


def compute_slab_groups(grid, prefix):
    gr = ""
    for i in grid.cell_data:
        if i.startswith(prefix):
            eids = np.where(grid.cell_data[i] > 0)[0] + 1
            gr += format_eset(i, eids)
    return gr


def nodebuffer(grid):
    nodes = np.column_stack((np.arange(1, len(grid.points) + 1), grid.points))
    logger.info(f"exporting {len(nodes)} nodes")
    return (
        "*node,nset=nall\n"
        + "\n".join([f"{int(n)},{x:f},{y:f},{z:f}" for n, x, y, z in nodes])
        + "\n"
    )


def element_buffer(grid):
    conn = grid.cells_dict
    logger.info(f"{conn}")
    extypes = [23]

    vtk_ccx = {23: "s8r"}

    buf = ""
    for tp in extypes:
        ccxtype = vtk_ccx[tp]
        for n, i in enumerate(conn[tp]):
            conn[tp][n] = np.array(i) + 1
            buf += f"*element,type={ccxtype},elset=e{n + 1}\n"
            buf += f"{n + 1},{','.join(map(str, conn[tp][n]))}\n"

    return buf


def orientation_buffer(grid, add_centers=False):
    shells = grid.extract_cells(grid.cells_dict.get(23, []))
    buf = ""
    x_dirs = shells.cell_data["x_dir"]
    y_dirs = shells.cell_data["y_dir"]
    centers = shells.cell_data["centers"]
    num_cells = shells.GetNumberOfCells()

    for n in range(num_cells):
        xdir, ydir = x_dirs[n], y_dirs[n]
        center = centers[n]
        buf += "*orientation,name=or%i,system=rectangular\n" % (n + 1)
        if add_centers:
            coords = np.concatenate([xdir + center, ydir + center, center], axis=0)
            buf += ",".join(format(k, ".4g") for k in coords.tolist()) + "\n3,0\n"
        else:
            coords = np.concatenate([xdir, ydir], axis=0)
            buf += ",".join(format(k, ".4g") for k in coords.tolist()) + "\n"

    return buf


def get_loadcases(mesh, multiplier=1.0, buckling=False):
    loadcases = {}

    for i in mesh.point_data:
        if i.startswith("lc_"):
            logger.info(f"loadcase {i}")
            multiplier = 1.0  # TODO fix for quadratic meshes
            if buckling:
                lbuf = f"** {i}\n*step\n*buckle\n5\n*cload\n"
            else:
                lbuf = f"** {i}\n*step\n*static\n*cload\n"

            ld = mesh.point_data[i] * multiplier
            for n, j in enumerate(ld):
                if j[0] ** 2 > 1e-8:
                    lbuf += "%i,1,%f\n" % (n + 1, j[0])
                if j[1] ** 2 > 1e-8:
                    lbuf += "%i,2,%f\n" % (n + 1, j[1])

            lbuf += "*node output,output=3d\nU\n*element output\nE,S\n*node print,nset=root,totals=yes\nrf\n*end step\n"

            loadcases[i] = lbuf

    return loadcases


def root_clamp(mesh):
    root = np.where(mesh.points[:, 2] == mesh.points[:, 2].min())
    lbuf = "*nset, nset=root\n"
    for n, j in enumerate(root[0]):
        lbuf += "%i" % (j + 1)
        if n % 16 == 0:
            lbuf += "\n"
        else:
            lbuf += ","
    lbuf += "\n"
    lbuf += "*boundary,op=new\nroot,1,3\n"
    return lbuf


def mesh2ccx(
    vtu,
    out="test.inp",
    matmap="temp/material_map.json",
    merge_adjacent_layers=True,
    zeroangle=False,
    single_step=False,
    quadratic=True,
    add_centers=False,
    force_isotropic=False,
    export_hyperworks=False,
    export_plygroups=False,
    buckling=False,
    meshonly=False,
    bondline=False,  # Added to accept bondline argument
):
    logger.info(f"** Running {vtu}")
    grid = pv.read(vtu)
    gr = grid.threshold(value=(1e-6, 1e9), scalars="thickness")
    gr.cell_data["centers"] = gr.cell_centers().points

    logger.info(f"** Exporting {gr.GetNumberOfCells()} elements")
    if quadratic:
        lf = vtk.vtkLinearToQuadraticCellsFilter()
        lf.SetInputData(gr)
        lf.Update()
        quad = lf.GetOutput()
        mesh = zero_midside_loads(pv.UnstructuredGrid(quad))
        mesh.save(vtu.replace(".vtu", "_quad.vtu"))
    else:
        mesh = pv.UnstructuredGrid(gr)

    buf = nodebuffer(mesh)

    if quadratic:
        buf += element_buffer(mesh)
    else:
        buf += element_buffer(mesh)

    if meshonly:
        of = out.replace(".inp", "_meshonly.inp")
        open(of, "w").write(buf)
        logger.info(f"written to {of}")
        return [of]

    buf += "*elset,elset=Eall,GENERATE\n%i,%i\n" % (1, mesh.GetNumberOfCells())

    if export_plygroups:
        plygroups, df = compute_ply_groups(mesh, "ply_")
        slabgroups = compute_slab_groups(mesh, "slab_thickness_")
        buf += plygroups + slabgroups

    plykeys = [i for i in mesh.cell_data if i.startswith("ply_")]
    plydat = np.stack([mesh.cell_data[i] for i in plykeys])
    materials = np.unique(plydat[:, :, 0])

    buf += orientation_buffer(mesh, add_centers)

    logger.info("made orientation buffer")
    matblock = material_db_to_ccx(materials, matmap=matmap, force_iso=force_isotropic)

    buf += "** START MATERIALS\n"
    buf += matblock
    buf += "** END MATERIALS\n"

    tic = time.perf_counter()

    blx = [
        make_shell_section(i, plydat[:, i, :], merge_adjacent_layers, zeroangle)
        for i in range(plydat.shape[1])
    ]
    logger.info("made shell sections")

    nplies = np.array([i[0] for i in blx])
    nplmax = nplies.max()
    npxid = np.where(nplies == nplmax)[0]
    logger.info(f"max number of plies: {nplmax}")
    # logger.info(f"associated stack \n{blx[npxid[0]][1]}")

    toc = time.perf_counter()
    logger.info(f"time spent creating shell sections {toc - tic:f}")
    logger.info(f"mesh {mesh}")
    comps = "".join(
        f"*shell section,composite,elset=e{n + 1},offset=-.5"
        + (f",orientation=or{n + 1}\n" if zeroangle else "\n")
        + i[1]
        for n, i in enumerate(blx)
    )
    buf += "** START SHELL SECTIONS\n"
    buf += comps
    buf += "** END SHELL SECTIONS\n"
    buf += root_clamp(mesh)

    loadcases = get_loadcases(mesh, buckling=buckling)

    logger.info(f"written loadcases {loadcases.keys()}")
    output_files = []
    if single_step:
        output = buf + "".join(loadcases.values())
        open(out, "w").write(output)
        logger.info(f"written ccx input file with all loadcases to {out}")
        output_files.append(out)
    else:
        for i in loadcases:
            output = buf + loadcases[i]
            of = out.replace(".inp", f"_{i}.inp")
            open(of, "w").write(output)
            output_files.append(of)
            logger.info(f"written ccx input file to {of}")

    if export_hyperworks:
        otb = out.replace(".inp", ".csv")
        df.to_csv(otb, index=False)
        logger.info(f"written plybook to hyperworks table {otb}")
        output_files.append(otb)

    return output_files
