#! /usr/bin/env python3

import pandas as pd
import pyvista
from io import StringIO
import multiprocessing
import numpy as np
import vtk
import argparse
import time


def findall(p, s):
    out = []
    i = s.find(p)
    out.append(i)
    while i != -1:
        i = s.find(p, i + 1)
        out.append(i)
    return out


def load(inp):
    if inp.find("CalculiX") != -1:  # read the nodes
        out = pd.read_fwf(
            StringIO(inp[inp.find("\n -1") :]),
            skiprows=1,
            names=["n", "x", "y", "z"],
            widths=[3, 10, 12, 12, 12],
            index_col=1,
        )
        return ("nodes", out.reset_index())
    elif inp.find("3C") > 1:  # elements
        e1 = pd.read_csv(
            StringIO(inp[inp.find("\n -1") :]),
            sep="\s+",
            skiprows=lambda x: (x % 3) != 0,
            names=range(11),
        )
        e2 = pd.read_csv(
            StringIO(inp[inp.find("\n -1") :]),
            sep="\s+",
            skiprows=lambda x: (x % 3) != 2,
            names=range(11),
        )
        return ("elements", pd.concat([e2.iloc[:, 1:], e1.iloc[:, 1:]], axis=1))
    elif inp.find("DISP") != -1:  # displacements
        out = pd.read_fwf(
            StringIO(inp[inp.find("\n -1") :]),
            skiprows=1,
            names=["n", "d1", "d2", "d3"],
            widths=[3, 10, 12, 12, 12],
            index_col=1,
        )
        return ("disp", out)
    elif inp.find("FORC") != -1:  # forces
        out = pd.read_fwf(
            StringIO(inp[inp.find("\n -1") :]),
            skiprows=1,
            names=["n", "f1", "f2", "f3"],
            widths=[3, 10, 12, 12, 12],
            index_col=1,
        )
        return ("force", out)
    elif inp.find("STRESS") != -1:  # stresses
        out = pd.read_fwf(
            StringIO(inp[inp.find("\n -1") :]),
            skiprows=1,
            names=["n", "sxx", "syy", "szz", "sxy", "syz", "szx"],
            widths=[3, 10, 12, 12, 12, 12, 12, 12],
            index_col=1,
        )
        return ("stress", out)
    elif inp.find("STRAIN") != -1:  # strains
        out = pd.read_fwf(
            StringIO(inp[inp.find("\n -1") :]),
            skiprows=1,
            names=["n", "exx", "eyy", "ezz", "exy", "eyz", "ezx"],
            widths=[3, 10, 12, 12, 12, 12, 12, 12],
            index_col=1,
        )
        return ("strain", out)

    return ("other", None)


def main():
    p = argparse.ArgumentParser(
        description="Fast translate of frd into vtu file using multiprocessing and pandas parsers, only tested on hex20"
    )
    p.add_argument("frd", help="Input file")
    p.add_argument("--output", default="__temp.vtu", help="Output file name")

    tic = time.perf_counter()
    args = p.parse_args()
    x = open(args.frd, "r").read()
    # find all block endings
    min3 = [0] + findall("\n -3", x)

    # pass the blocks to a process each
    p = multiprocessing.Pool()
    o = dict(p.map(load, [x[i[0] : i[1]] for i in zip(min3, min3[1:])]))

    # create a map between ccx node id and vtk
    idmap = pd.DataFrame(o["nodes"].index, o["nodes"]["index"])

    # create new connectivity table
    conn = idmap.loc[o["elements"].values.flatten()].values
    conn.resize(o["elements"].shape)

    # renumber the hex20 element nodes according to vtk definition
    conn = np.concatenate([conn[:, :12], conn[:, 16:], conn[:, 12:16]], axis=1)

    # use the pyvista dict mesh build interface
    ogrid = pyvista.UnstructuredGrid(
        {vtk.VTK_QUADRATIC_HEXAHEDRON: conn}, o["nodes"][["x", "y", "z"]].values
    )
    for i in ["disp", "force", "stress", "strain"]:
        ogrid.point_data[i] = o[i].values[:, 1:]

    # add mises strain and stress
    for s in [("stress", "s"), ("strain", "e")]:
        ss = o[s[0]]
        xx, yy, zz, xy, yz, zx = [
            ss["%s%s" % (s[1], i)] for i in ["xx", "yy", "zz", "xy", "yz", "zx"]
        ]
        ogrid.point_data["mises_%s" % s[0]] = (1.0 / 2.0 ** 0.5) * (
            (xx - yy) ** 2
            + (yy - zz) ** 2
            + (zz - xx) ** 2
            + 6.0 * (xy ** 2 + yz ** 2 + zx ** 2)
        ) ** 0.5

    ogrid.save(args.output)

    print("written to ", args.output)
    toc = time.perf_counter()
    print("time for translation %fs" % (toc - tic))


if __name__ == "__main__":
    main()
