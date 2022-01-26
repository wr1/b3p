#! /usr/bin/env python3

import pyvista
import numpy as np
import multiprocessing
import vtk
import argparse
from functools import partial
import time
import json
import os
import yaml


def make_shell_section(plyarray):
    plies = []
    # if material in next ply is the same, add it to the previous ply,
    # NOTE this does not take ply angle into account so it only works
    # for multidirectional fibre mats at angle==0 for now
    for j in plyarray[np.where(plyarray[:, 1] > 1e-6)]:
        # ply_dat = g.cell_data[j][n, :]
        if len(plies) > 0 and plies[-1][1] == j[0]:
            plies[-1][0] += j[1] * 1e-3
        else:
            plies.append([j[1] * 1e-3, j[0]])

    comp = ""
    for i in plies:
        comp += "%f,,m%i,or1\n" % tuple(i)
    return comp


def material_db_to_ccx(grid, materials):
    gdir = os.path.dirname(grid)
    mm_name = os.path.join(gdir, "material_map.json")
    # 2 files are relevant, a material map that maps the material ID (integer)
    # in the VTK file to a key (string) in the material database and the material
    # database itself
    mat_db = None
    if os.path.isfile(mm_name):  # check if the material map file is there
        mm = json.load(open(mm_name, "r"))
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


def main():
    p = argparse.ArgumentParser(description="Translate a blade vtk mesh into ccx")
    p.add_argument(
        "grid",
        help="Grid file (vtu) including shear web(s) and shell (assuming you want to simulate a blade)",
    )
    p.add_argument("--out", default="test.inp", help="Output file name (input for ccx)")
    args = p.parse_args()

    grid = args.grid

    re = vtk.vtkXMLUnstructuredGridReader()
    re.SetFileName(grid)
    re.Update()
    gr = re.GetOutput()

    lf = vtk.vtkLinearToQuadraticCellsFilter()
    lf.SetInputData(re.GetOutput())
    lf.Update()

    quad = lf.GetOutput()

    g = pyvista.UnstructuredGrid(quad)

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

    # write orientation TODO match with element orientation, for now just align with z-axis
    for n, i in enumerate(conn):
        buf += "*orientation,name=or%i,system=rectangular\n" % (n + 1)
        buf += "0,0,1,0,1,0\n"
        break

    plydat = np.stack(g.cell_data[i] for i in g.cell_data if i.startswith("ply_"))
    # get all materials of all plies
    materials = np.unique(plydat[:, :, 0])

    matblock = material_db_to_ccx(grid, materials)

    buf += matblock

    tic = time.perf_counter()
    p = multiprocessing.Pool()

    blx = p.map(make_shell_section, [plydat[:, i, :] for i in range(plydat.shape[1])])
    toc = time.perf_counter()
    print("time spent creating shell sections %f" % (toc - tic))

    # unique_laminates = np.unique(blx)

    comps = ""
    for n, i in enumerate(blx):
        comps += "*shell section, composite, elset=e%i,offset=-.5\n%s" % (n + 1, i)

    buf += comps

    ## this minimizes the number of different composite sections, however this is not a big factor,
    ## what causes the model size explosion is the number of layers in the composite cards
    # elset = ""
    # comps = ""
    # for n, i in enumerate(unique_laminates):
    #     occ = np.where(np.array(blx) == i)[0]
    #     elset += "*elset,elset=e%i\n" % (n + 1)
    #     for m, j in enumerate(occ):
    #         elset += "%i" % (j + 1)
    #         if m % 16 == 15 or j == occ[-1]:
    #             elset += "\n"
    #         else:
    #             elset += ","
    #     comps += "*shell section, composite, elset=e%i,offset=0\n%s" % (n + 1, i)
    # buf += elset + comps

    buf += "*step\n*static\n"

    mid = np.where(
        (g.points[:, 2] > g.points[:, 2].max() * 0.7)
        & (g.points[:, 2] < g.points[:, 2].max() * 0.8)
    )

    loadcases = {}

    for i in g.point_data:
        if i.startswith("lc_"):
            lbuf = ""
            ld = g.point_data[i]
            for n, j in enumerate(ld):
                if j[0] ** 2 > 1e-8:
                    lbuf += "%i,1,%f\n" % (n + 1, j[0])
                if j[1] ** 2 > 1e-8:
                    lbuf += "%i,2,%f\n" % (n + 1, j[1])

            loadcases[i] = "*cload\n" + lbuf

    for i in loadcases:
        buf += loadcases[i]
        break  # only add the first loadcase

    frc = "*cload\n"
    for i in mid[0]:
        frc += "%i,1,1e6\n" % (i + 1)

    root = np.where(g.points[:, 2] == g.points[:, 2].min())
    bcs = "*boundary,op=new\n"
    for i in root[0]:
        bcs += "%i,1,3\n" % (i + 1)

    buf += bcs + frc + "*node file,output=3d\nU,RF\n*EL FILE\nS,E\n" + "*end step\n"

    open(args.out, "w").write(buf)
    print("written ccx input file to %s" % args.out)


if __name__ == "__main__":
    main()