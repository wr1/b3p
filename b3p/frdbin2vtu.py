#!/usr/bin/env python
import numpy as np
import re
import pandas as pd
import pyvista as pv
import vtk
import fire
import time


def split_blocks(buf):
    patterns = [b" 2C", b"  3C", b"PSTEP", b"DISP", b"STRAIN", b"FORC", b" 9999"]
    locations = {
        pattern.decode("ascii").strip(): [m.start() for m in re.finditer(pattern, buf)]
        for pattern in patterns
    }

    return locations


def is_ascii(byte):
    return 0x20 <= byte <= 0x7E or byte in {0x0A, 0x0D}


def extract_binary_data(block):
    start, lastn = 0, 0
    end = len(block)
    # Skip the initial ASCII part
    while start < end and is_ascii(block[start]):
        if chr(block[start]) == "\n":
            lastn = start
        start += 1

    # Skip the trailing ASCII part
    trailing_start = end - 1
    lastspace = end - 1
    while trailing_start > start and is_ascii(block[trailing_start]):
        if chr(block[trailing_start]) == " ":
            lastspace = trailing_start
        trailing_start -= 1

    # Return the binary data in the middle
    return block[lastn + 1 : lastspace]


def load_block(buf, start, end, dtype):
    # bin = buf[start:end]
    # print(bin[:400], bin[-20:])
    bind = extract_binary_data(buf[start:end])

    # print(bind[:20], bind[-20:])
    # print(f"{bind[0]}")
    return pd.DataFrame(np.frombuffer(bind, dtype=dtype))


def frdbin2vtu(file_path):
    starttime = time.time()
    print(f"Converting {file_path}")
    buf = open(file_path, "rb").read()
    lcs = split_blocks(buf)

    nodes = load_block(
        buf,
        lcs["2C"][0],
        lcs["3C"][0],
        np.dtype([("i", "i4"), ("x", "f8"), ("y", "f8"), ("z", "f8")]),
    )

    # print(nodes[0])
    # exit()
    conn = load_block(
        buf,
        lcs["3C"][0],
        lcs["PSTEP"][0],
        np.dtype([("i%i" % i, "i4") for i in range(24)]),
    )

    disp = load_block(
        buf,
        lcs["PSTEP"][0],
        lcs["PSTEP"][1],
        np.dtype([("i", "i4"), ("x", "f4"), ("y", "f4"), ("z", "f4")]),
    )
    stype = np.dtype(
        [
            ("i", "i4"),
            ("xx", "f4"),
            ("yy", "f4"),
            ("zz", "f4"),
            ("xy", "f4"),
            ("yz", "f4"),
            ("xz", "f4"),
        ]
    )
    stress = load_block(
        buf,
        lcs["PSTEP"][1],
        lcs["PSTEP"][2],
        # 0,
        # 2,
        stype,
    )

    strain = load_block(
        buf,
        lcs["PSTEP"][2],
        lcs["PSTEP"][3] if len(lcs["PSTEP"]) > 3 else lcs["9999"][0],
        stype,
    )

    # reorder the connectivity array, and map it from ccx node id to vtk node id
    nz = np.zeros(nodes["i"].max() + 1, dtype=int)
    nz[nodes["i"]] = nodes.index

    hx = conn.values[:, 4:]
    hx = np.hstack([hx[:, :12], hx[:, 16:], hx[:, 12:16]])
    hshape = hx.shape

    hx = nz[hx.flatten()].reshape(hshape)

    # use the pyvista dict mesh build interface
    ogrid = pv.UnstructuredGrid(
        {vtk.VTK_QUADRATIC_HEXAHEDRON: hx},
        nodes[["x", "y", "z"]].values,
    )

    # add cell data
    ogrid.cell_data["material"] = conn["i3"].values
    ogrid.cell_data["ccx_id"] = conn["i0"].values
    ogrid.cell_data["i2"] = conn["i2"].values
    ogrid.cell_data["i1"] = conn["i1"].values

    # add point data
    tm = "1.00000"
    ogrid.point_data[f"disp_{tm}"] = disp[["x", "y", "z"]].values
    ogrid.point_data[f"strain_{tm}"] = strain[
        ["xx", "yy", "zz", "xy", "yz", "xz"]
    ].values
    ogrid.point_data[f"stress_{tm}"] = stress[
        ["xx", "yy", "zz", "xy", "yz", "xz"]
    ].values
    ogrid.point_data["ccx_id"] = nodes["i"].values
    of = file_path.replace(".frd", ".vtu")
    ogrid.save(of)
    print(f"Saved {of}")
    endtime = time.time()
    print(f"Elapsed time: {endtime - starttime} seconds")


def main():
    fire.Fire(frdbin2vtu)


if __name__ == "__main__":
    main()
