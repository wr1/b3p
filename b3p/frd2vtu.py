#! /usr/bin/env python3

import pandas as pd
import pyvista
from io import StringIO
import multiprocessing
import numpy as np
import vtk
import fire
import time

import regex


def findall(p, s):
    return [m.start() for m in regex.finditer(p, s, overlapped=True)]


def load_frd(inp):
    if inp.find("PSTEP") != -1:
        start_step_line = inp.find("100CL")
        end_step_line = inp.find("\n", start_step_line)
        _, solution_type, timestep = inp[start_step_line:end_step_line].split()[:3]

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
        return (
            "elements",
            pd.concat([e2.iloc[:, 1:], e1.iloc[:, 1:]], axis=1),
        )
    elif inp.find("DISP") != -1:  # displacements
        out = pd.read_fwf(
            StringIO(inp[inp.find("\n -1") :]),
            skiprows=1,
            names=["n", "d1", "d2", "d3"],
            widths=[3, 10, 12, 12, 12],
            index_col=1,
        ).astype("float32")
        return (f"disp_{timestep}", out)
    elif inp.find("FORC") != -1:  # forces
        out = pd.read_fwf(
            StringIO(inp[inp.find("\n -1") :]),
            skiprows=1,
            names=["n", "f1", "f2", "f3"],
            widths=[3, 10, 12, 12, 12],
            index_col=1,
        ).astype("float32")
        return (f"force_{timestep}", out)
    elif inp.find("STRESS") != -1:  # stresses
        out = pd.read_fwf(
            StringIO(inp[inp.find("\n -1") :]),
            skiprows=1,
            names=["n", "sxx", "syy", "szz", "sxy", "syz", "szx"],
            widths=[3, 10, 12, 12, 12, 12, 12, 12],
            index_col=1,
        ).astype("float32")
        return (f"stress_{timestep}", out)
    elif inp.find("STRAIN") != -1:  # strains
        out = pd.read_fwf(
            StringIO(inp[inp.find("\n -1") :]),
            skiprows=1,
            names=["n", "exx", "eyy", "ezz", "exy", "eyz", "ezx"],
            widths=[3, 10, 12, 12, 12, 12, 12, 12],
            index_col=1,
        ).astype("float32")
        return (f"strain_{timestep}", out)

    return ("other", None)


def frd2vtu(frd, output=None, multi=False):
    if not output:
        output = frd.replace(".frd", ".vtu")

    print("running frd2vtu on file: ", frd)
    tic = time.perf_counter()
    x = open(frd, "r").read()
    # find all block endings
    min3 = [0] + findall("\n -3", x)

    # pass the blocks to a process each
    if multi:
        p = multiprocessing.Pool()
        o = p.map(load_frd, [x[i[0] : i[1]] for i in zip(min3, min3[1:])])
    else:
        o = [load_frd(x[i[0] : i[1]]) for i in zip(min3, min3[1:])]

    o = dict(o)

    # create a map between ccx node id and vtk
    idmap = pd.DataFrame(o["nodes"].index, o["nodes"]["index"])

    # create new connectivity table
    el = o["elements"]
    # hex elements
    nna = el.isna().sum(axis=1).values
    hx = el.loc[nna == 0, :].values.astype(int)
    hx = np.hstack([hx[:, :12], hx[:, 16:], hx[:, 12:16]])

    # wedge elements
    wd = el.loc[nna == 5, :].values[:, :-5].astype(int)
    wd = np.hstack([wd[:, :9], wd[:, 12:], wd[:, 9:12]])

    hshape = hx.shape
    hx = idmap.loc[hx.flatten()].values.reshape(hshape)

    wshape = wd.shape
    wd = idmap.loc[wd.flatten()].values.reshape(wshape)

    # use the pyvista dict mesh build interface
    ogrid = pyvista.UnstructuredGrid(
        {vtk.VTK_QUADRATIC_HEXAHEDRON: hx, vtk.VTK_QUADRATIC_WEDGE: wd},
        o["nodes"][["x", "y", "z"]].values,
    )

    for i in o:
        if any(
            i.startswith(prefix) for prefix in ["disp", "force", "stress", "strain"]
        ):
            ogrid.point_data[i] = o[i].values[:, 1:]

    # add mises strain and stress

    # for s in [("stress", "s"), ("strain", "e")]:
    #     ss = o[s[0]]
    #     xx, yy, zz, xy, yz, zx = [
    #         ss[f"{s[1]}{i}"] for i in ["xx", "yy", "zz", "xy", "yz", "zx"]
    #     ]
    #     ogrid.point_data[f"mises_{s[0]}"] = (1.0 / 2.0**0.5) * (
    #         (xx - yy) ** 2
    #         + (yy - zz) ** 2
    #         + (zz - xx) ** 2
    #         + 6.0 * (xy**2 + yz**2 + zx**2)
    #     ) ** 0.5

    ogrid.save(output)

    print(f"** written to {output}")
    toc = time.perf_counter()
    print("time for translation %fs" % (toc - tic))
    return ogrid


if __name__ == "__main__":
    fire.Fire(frd2vtu)
