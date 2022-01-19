#! /usr/bin/env python3

import pyvista
import numpy as np
import multiprocessing
import vtk
import argparse
from functools import partial
import time


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
        comp += "%f,,mat%i,or1\n" % tuple(i)
    return comp


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

    # TODO proper material definitions defined in yaml file...
    mod = {1: 100e9, 2: 40e9, 3: 10e9, 4: 20e9, 5: 0.2e9, 6: 0.2e9, 7: 10e9, 8: 0.2e9}
    matblock = ""
    for i in materials:
        if i > 1e-6:
            matblock += "*material,name=mat%i\n*elastic,type=iso\n" % i
            matblock += "%f,.3\n" % mod[i]

    buf += matblock

    tic = time.perf_counter()
    p = multiprocessing.Pool()

    blx = p.map(make_shell_section, [plydat[:, i, :] for i in range(plydat.shape[1])])
    toc = time.perf_counter()
    print("time spent creating shell sections %f" % (toc - tic))

    # unique_laminates = np.unique(blx)

    comps = ""
    for n, i in enumerate(blx):
        comps += "*shell section, composite, elset=e%i,offset=-1\n%s" % (n + 1, i)

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


if __name__ == "__main__":
    main()