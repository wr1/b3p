#! /usr/bin/env python3

import pyvista as pv
import numpy as np
import multiprocessing
import vtk
import argparse
from functools import partial
import time
import json
import os
import yaml


def make_shell_section(inp):
    elem_id, plyarray = inp
    plies = []
    # if material in next ply is the same, add it to the previous ply,
    # NOTE this does not take ply angle into account so it only works
    # for multidirectional fibre mats at angle==0 for now
    for j in plyarray[np.where(plyarray[:, 1] > 1e-6)]:
        # ply_dat = g.cell_data[j][n, :]
        if plies and plies[-1][1] == j[0]:
            plies[-1][0] += j[1] * 1e-3
        else:
            plies.append([j[1] * 1e-3, j[0]])

    return "".join("%f,,m%i,or%i\n" % tuple(i + [elem_id + 1]) for i in plies)


def material_db_to_ccx(grid, materials, matmap=None):
    "find the material db and write properties to a ccx block"

    # 2 files are relevant, a material map that maps the material ID (integer)
    # in the VTK file to a key (string) in the material database and the material
    # database itself
    mat_db = None
    if os.path.isfile(matmap):  # check if the material map file is there
        gdir = os.path.dirname(matmap)
        mm = json.load(open(matmap, "r"))
        if "matdb" in mm:  # check if the material map file points to a material db
            mat_db = yaml.load(
                open(os.path.join(gdir, mm["matdb"])), Loader=yaml.CLoader
            )
        else:
            exit(
                "material map available, but no link to material db, need matdb definition to do FEA"
            )
    else:
        exit("no material map defined")

    mm_inv = {v: k for k, v in mm.items()}
    matblock = ""
    for i in materials:
        if i > 1e-6:
            material_properties = mat_db[mm_inv[int(i)]]

            if "vf" in material_properties:
                print(material_properties["name"], "is assumed to be orthotropic")
                C = np.array(material_properties["C"])
                # https://github.com/rsmith-nl/lamprop/blob/410ebfef2e14d7cc2988489ca2c31103056da38f/lp/text.py#L96
                # https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node193.html
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
                D *= 1e6
                matblock += (
                    f"{D[0,0]:.4g},{D[0,1]:.4g},{D[1,1]:.4g},"
                    + f"{D[0,2]:.4g},{D[1,2]:.4g},{D[2,2]:.4g},"
                    + f"{D[3,3]:.4g},{D[4,4]:.4g},\n"
                    + f"{D[5,5]:.4g},293\n"
                )
            else:
                print(material_properties["name"], "is assumed to be isotropic")
                nu = min(0.45, max(0.1, float(material_properties["nu"])))  # min()
                E = 1e6 * float(
                    material_properties["Ex"]
                    if "Ex" in material_properties
                    else material_properties["E"]
                )
                matblock += "*material,name=m%i\n*elastic,type=iso\n" % i
                matblock += f"{E:.4g},{nu:.4g},293\n"

    return matblock


def format_eset(name, eids, prefix):
    out = f"*elset,elset={name.replace(prefix,'')}\n"
    for i in range(len(eids)):
        out += f"{eids[i]}"
        out += "\n" if (i % 16 == 15) else ","
    if out[-1] == ",":
        out = out[:-1] + "\n"
    return out


def compute_slab_groups(grid):
    gr = ""
    prefix = "slab_thickness_"
    for i in grid.cell_data:
        if i.startswith(prefix):
            eids = np.where(grid.cell_data[i] > 0)[0] + 1
            gr += format_eset(i, eids, prefix)

    return gr


def main():
    p = argparse.ArgumentParser(description="Translate a blade vtk mesh into ccx")
    p.add_argument(
        "grid",
        help="Grid file (vtu) including shear web(s) and shell (assuming you want to simulate a blade)",
    )
    p.add_argument("--out", default="test.inp", help="Output file name (input for ccx)")
    p.add_argument("--matmap", default="temp/material_map.json")
    args = p.parse_args()

    grid = args.grid

    gr = pv.read(grid)
    glin = pv.UnstructuredGrid(gr)

    lf = vtk.vtkLinearToQuadraticCellsFilter()
    lf.SetInputData(gr)  # re.GetOutput())
    lf.Update()
    quad = lf.GetOutput()
    g = pv.UnstructuredGrid(quad)

    # export the nodes
    buf = "*node,nset=nall\n"
    for n, i in enumerate(g.points):
        buf += "%i,%f,%f,%f\n" % tuple([n + 1] + list(i))

    conn = g.cell_connectivity.reshape(
        (g.GetNumberOfCells(), int(g.cell_connectivity.shape[0] / g.GetNumberOfCells()))
    )

    # export the elements
    for n, i in enumerate(conn):
        buf += "*element,type=s8r,elset=e%i\n" % (n + 1)
        buf += "%i,%i,%i,%i,%i,%i,%i,%i,%i\n" % tuple([n + 1] + list(i + 1))

    buf += "*elset,elset=Eall,GENERATE\n%i,%i\n" % (1, n + 1)

    buf += compute_slab_groups(g)

    # write orientation TODO match with element orientation, for now just align with z-axis
    for n, i in enumerate(conn):
        xdir, ydir = g.cell_data["x_dir"][n], g.cell_data["y_dir"][n]
        buf += "*orientation,name=or%i,system=rectangular\n" % (n + 1)
        buf += "%.4g,%.4g,%.4g,%.4g,%.4g,%.4g\n" % tuple(xdir.tolist() + ydir.tolist())

    plydat = np.stack(g.cell_data[i] for i in g.cell_data if i.startswith("ply_"))
    # get all materials of all plies
    materials = np.unique(plydat[:, :, 0])

    matblock = material_db_to_ccx(grid, materials, matmap=args.matmap)

    buf += matblock

    tic = time.perf_counter()
    p = multiprocessing.Pool()

    blx = p.map(
        make_shell_section, [(i, plydat[:, i, :]) for i in range(plydat.shape[1])]
    )
    toc = time.perf_counter()
    print("time spent creating shell sections %f" % (toc - tic))

    comps = "".join(
        "*shell section, composite, elset=e%i,offset=-.5\n%s" % (n + 1, i)
        for n, i in enumerate(blx)
    )
    buf += comps

    buf += "*step\n*static\n"

    # mid = np.where(
    #     (g.points[:, 2] > g.points[:, 2].max() * 0.7)
    #     & (g.points[:, 2] < g.points[:, 2].max() * 0.8)
    # )

    loadcases = {}

    for i in g.point_data:
        if i.startswith("lc_"):
            # forces are interpolated to midside nodes, causing the sum of forces to be off,
            # compute a multiplier from the sum of the forces in the linear model here
            multiplier = glin.point_data[i].sum() / g.point_data[i].sum()
            lbuf = ""
            ld = g.point_data[i] * multiplier
            for n, j in enumerate(ld):
                if j[0] ** 2 > 1e-8:
                    lbuf += "%i,1,%f\n" % (n + 1, j[0])
                if j[1] ** 2 > 1e-8:
                    lbuf += "%i,2,%f\n" % (n + 1, j[1])

            loadcases[i] = "*cload\n" + lbuf

    root = np.where(g.points[:, 2] == g.points[:, 2].min())
    bcs = "*boundary,op=new\n"
    for i in root[0]:
        bcs += "%i,1,3\n" % (i + 1)

    endfile = (
        bcs
        + "*node file,output=3d\nU,RF\n*EL FILE\nS,E\n*node print,nset=nall\nrf\n"
        + "*end step\n"
    )

    # write a full ccx file for each loadcase, assuming parallel execution
    for i in loadcases:
        print(i)
        output = buf + loadcases[i] + endfile
        of = args.out.replace(".inp", f"_{i}.inp")
        open(of, "w").write(output)
        print(f"written ccx input file to {of}")


if __name__ == "__main__":
    main()
