#! /usr/bin/env python3

import numpy as np

from ruamel.yaml import YAML
from itertools import chain, zip_longest
import pickle
import json
import copy as cp
from sympy import symbols, sympify
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def plyify(r, t, ply_thickness, reverse=False):
    """Fill thickness distribution using plies.

    params:
        r (list): radius distribution
        t (list): thickness distribution
        ply_thickness (float): thickness of the plies
        reverse (bool): flag determining whether the plies are numbered from the
            top or bottom of the stack
        returns:    list of plies in the form [rmin, rmax]"""
    active = []
    done = []
    np = len(r)
    for i in range(np):
        while t[i] > len(active) * ply_thickness:
            active.append([r[i], -1])
        while t[i] <= (len(active) - 1) * ply_thickness:
            pop = active.pop()
            pop[-1] = r[i]
            done.append(pop)
    for j in active:
        j[-1] = r[i]
        done.append(j)
    for i in done:
        i[0] = round(i[0] * 1e3, -1) * 1e-3
        i[1] = round(i[1] * 1e3, -1) * 1e-3
    # reverse the stack so that longest plies get the lowest ply numbers
    return done if reverse else list(reversed(done))


def expand_chamfered_core(core):
    rr, tr = np.array(core["slab"]).T

    coordinates = []
    added_cores = []
    for chmf in core["chamfers"]:
        ns = chmf["nstep"] if "nstep" in chmf else 4

        idx = chmf["id"]

        fr = np.array(chmf["ratio"])

        intrr = np.sort(np.unique(np.hstack([rr, fr[:, 0]])))

        thickness = np.interp(intrr, rr, tr)

        ratio = np.interp(intrr, fr[:, 0], fr[:, 1])

        distance = 1e-3 * ratio * thickness

        # chamfer_step = "chstep"
        added_coordinates = {}
        for i in range(0, ns + 1):
            key = f"{idx[0]}_c{i}"
            added_coordinates[key] = {
                "base": idx[0],
                "points": np.stack(
                    [intrr, (-1.0 if idx[1] == 0 else 1.0) * i * distance / ns]
                ).T.tolist(),
            }

        coordinates.append(added_coordinates)

        for j in range(0, ns):
            thic = thickness * ((j + 1.0) / (ns + 1.0))
            nc = cp.deepcopy(core)
            del nc["chamfers"]

            nc["slab"] = np.stack([intrr, thic]).T.tolist()

            bc = nc["cover"][f"{idx[0]}"]
            if idx[1] == 0:
                nc["cover"][f"{idx[0]}_c{j}"] = bc
                nc["cover"][f"{idx[0]}_c{j + 1}"] = [-1, bc[0], 0]

            # chamfers from 2nd coordinate
            elif idx[1] == 1:
                nc["cover"][f"{idx[0]}_c{j}"] = bc
                nc["cover"][f"{idx[0]}_c{j + 1}"] = [bc[1], -bc[0], 0]

            del nc["cover"][f"{idx[0]}"]

            added_cores.append(nc)

        core["cover"][key] = cp.deepcopy(core["cover"][idx[0]])

        if idx[0] in core["cover"]:
            del core["cover"][idx[0]]

    del core["chamfers"]
    return coordinates, added_cores


def expand_chamfered_cores(bldict, outputfile: Path = None):
    """ "Expand chamfered cores into laminate strips and write out to a new YAML file."""
    dctcopy = cp.deepcopy(bldict)
    for i in bldict["laminates"]["slabs"]:
        slab = dctcopy["laminates"]["slabs"][i]
        if "chamfers" in slab:
            coordinates, cores = expand_chamfered_core(slab)
            for j in range(len(cores)):
                slabname = i + f"_ch{j}"
                dctcopy["laminates"]["slabs"][slabname] = cores[j]

            for j in range(len(coordinates)):
                for k in coordinates[j]:
                    dctcopy["mesh"]["coordinates"][k] = coordinates[j][k]

    yaml = YAML()
    yaml.dump(dctcopy, outputfile.open("w"))
    if outputfile:
        logger.info(f"written expanded chamfered cores to {outputfile}")
    else:
        logger.warning(
            "No output file specified, not writing expanded chamfered cores."
        )
    return dctcopy


def ply_stack(r, t, t_ply=1.0, reverse=False, subdivisions=5000, material=11):
    """
    Taking in a r,t thickness over radius distribution, split it up in
    plies, note that it rounds down in thickness (if the last ply does not fit
    in the thickness it does not go on), so if you have thick plies,
    make sure that the mm thickness is slightly up from a whole number of plies.

    :param r: radius distribution
    :param t: thickness distribution
    :param t_ply: thickness of the plies
    :param reverse: flag determining whether the plies are numbered from the
        top or bottom of the stack
    :param subdivisions: number of subdivisions
    :param material: material number
    """
    x = np.linspace(min(r), max(r), subdivisions)
    y = np.interp(x, list(r), t)
    st = plyify(x, y, t_ply, reverse)
    return [[material, t_ply] + i for i in st]


def coreblock(r, t, subdivisions=200, material=11):
    """Make a thickness distribution into a block divisions.

    :param r: radius distribution
    :param t: thickness distribution
    :param subdivisions: number of subdivisions
    :param material: material number"""
    assert len(r) == len(t)
    lr = len(r)
    x = np.array(sorted(list(np.linspace(min(r), max(r), int(subdivisions))) + list(r)))
    np.interp(0.5 * (x[:-1] + x[1:]), list(r), t)

    stack = []
    for i in range(lr - 1):
        rmin, rmax = r[i], r[i + 1]
        tmin, tmax = t[i], t[i + 1]
        tt = tmin if tmin == tmax else 0.5 * (tmin + tmax)
        if tt > 0:
            stack.append([material, tt, rmin, rmax])

    return stack


def number_stack(stack, splitstack, key, increment):
    """Create ply numbering for the plies in the stack, depends on the start key
    and increment as well as the above/below ratio (splitstack).

    :param stack: list of plies
    :param splitstack: ratio of above/below plies
    :param key: start key from top and bottom
    :param increment: increment, number by which to increment the key (can be used to interleave plies)
    """
    sp = np.nan_to_num(splitstack.astype(float))
    if sp.sum() != 1.0:
        exit(f"sum {splitstack} not 1.")

    ky = np.nan_to_num(key.astype(float)).astype(int)
    inc = np.nan_to_num(increment.astype(float), nan=1).astype(int)

    ss = (len(stack) * sp).round().astype(int)
    st, sb = np.array(range(ss[0])), np.array(range(ss[1]))
    stt = ky[0] + st * inc[0]
    sbt = ky[1] + sb * inc[1]
    return [i for i in chain.from_iterable(zip_longest(stt, sbt)) if i is not None]


from sympy import lambdify


def get_coverage(slab, datums, rr):
    """
    Get the coverage of the slab and insert the datums.
    :param slab: slab dictionary
    :param datums: datums dictionary
    :param rr: radius distribution
    """
    assert "cover" in slab
    cov = slab["cover"]

    # Create symbolic variables for all datum keys
    symbols_dict = {k: symbols(k) for k in datums}

    for d in cov:
        cc = cov[d]
        # Parse the expressions symbolically
        idd = [sympify(cc[0], locals=symbols_dict), sympify(cc[1], locals=symbols_dict)]

        for j in range(2):
            expr = idd[j]  # Symbolic expression

            # Replace symbols with interpolated numerical arrays
            substitutions = {}
            for key, sym in symbols_dict.items():
                if sym in expr.free_symbols:  # Check if the symbol is present
                    xy = np.array(datums[key]["xy"])
                    dst = np.interp(
                        rr,
                        xy[:, 0] / datums[key]["scalex"],
                        xy[:, 1] * datums[key]["scaley"],
                    )
                    substitutions[sym] = dst  # Store the substitution separately

            # Perform the substitutions numerically after creating a callable function
            func = lambdify(list(substitutions.keys()), expr, modules="numpy")
            idd[j] = func(*substitutions.values())  # Pass the numerical arrays

        cov[d][0], cov[d][1] = idd  # Assign the evaluated results

    return cov


def add_bondline_material(matdb, material_map):
    """
    Add a bondline material to the material map, if it is not already there,
    this is needed because the 2d mesher adds elements with -1 material ID to connect the webs to the shell

    params:
        matdb (dict): material database
        material_map (dict): material map
        returns:
            material_map (dict): material map with bondline material added
    """
    possible_glue_names = ["adhesive", "bonding"]
    for i in matdb:
        for j in possible_glue_names:
            if matdb[i]["name"].find(j) != -1:
                matkey = i
                material_map["bondline"] = -1

    if "bondline" in material_map:
        matdb["bondline"] = cp.deepcopy(matdb[matkey])
        matdb["bondline"]["name"] = "bondline"

    return material_map


def export_matdb(blade, material_map, material_map_file: Path = None):
    """Export material database to JSON."""
    # workdir = blade["general"]["workdir"]
    # matmap = os.path.join(workdir, "drape", "material_map.json")
    matmap = material_map_file

    if "bondline" in blade["mesh"] and "material" in blade["mesh"]["bondline"]:
        maxid = max(v for v in material_map.values() if isinstance(v, int))
        material_map[blade["mesh"]["bondline"]["material"]] = maxid + 1

    material_map = add_bondline_material(blade["materials"], material_map)

    with open(matmap, "w") as f:
        json.dump({"matdb": blade["materials"], "map": material_map}, f)

    logger.info(f"written material map to {matmap}")


def export_plybook(stacks, outputfile):
    pickle.dump(stacks, open(outputfile, "wb"))
    logger.info(f"written to {outputfile}")


def lamplan2plies(blade, outputfile="__plybook.pck"):
    """Convert a lamplan to a list of plies."""
    root_radius = blade["planform"]["z"][0][1]

    tip_radius = blade["planform"]["z"][-1][1]

    slabs = blade["laminates"]["slabs"]

    datums = blade["laminates"]["datums"] if "datums" in blade["laminates"] else {}

    allstacks = []

    # use a multiple of lamplan length as radius grid to interpolate geometric variables to
    n_s = round(tip_radius * 4)

    radius = np.linspace(0, tip_radius - root_radius, n_s)
    radius_relative = np.linspace(0, 1, n_s)
    material_map = {}

    for i in slabs:
        material = slabs[i]["material"]
        if material not in material_map:
            material_map[material] = len(material_map) + 1

        grid = slabs[i]["grid"]

        cover = get_coverage(slabs[i], datums, radius_relative)

        draping = "plies" if "draping" not in slabs[i] else slabs[i]["draping"]

        ply_thickness = float(
            1.0 if "ply_thickness" not in slabs[i] else slabs[i]["ply_thickness"]
        )

        r, t = np.array(slabs[i]["slab"]).T

        r *= (
            tip_radius - root_radius if "rscale" not in slabs[i] else slabs[i]["rscale"]
        )
        t *= ply_thickness

        if draping == "blocks":
            stack = coreblock(r, t, material=material_map[material])
        elif draping == "plies":
            stack = ply_stack(r, t, ply_thickness, material=material_map[material])

        # assigns keys to the plies in the stack depending on how the stack is split up
        # what increment the ply keys are stacked with (non-1 increment allowing interleaving of plies)
        stack_numbering = number_stack(
            stack,
            (
                np.array([1, 0])
                if "splitstack" not in slabs[i]
                else np.array(slabs[i]["splitstack"])
            ),
            (
                np.array([0, 4000])
                if "key" not in slabs[i]
                else np.array(slabs[i]["key"])
            ),
            (
                np.array([1, -1])
                if "increment" not in slabs[i]
                else np.array(slabs[i]["increment"])
            ),
        )
        allstacks.append(
            {
                "name": i.strip(),
                "grid": grid.strip(),
                "cover": cover,
                "numbering": stack_numbering,
                "stack": stack,
                "r": radius,
            }
        )
    export_matdb(
        blade,
        material_map,
        material_map_file=Path(outputfile).parent / "material_map.json",
    )
    export_plybook(allstacks, outputfile=outputfile)
    return allstacks


def slab2plybook(yamlfile, outputfile="__lamplan.pck"):
    blade = yaml.safe_load(open(yamlfile, "r"))
    lamplan2plies(blade, outputfile)
    logger.info(f"written plydrape to {outputfile}")
