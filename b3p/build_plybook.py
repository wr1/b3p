#! /usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import yaml
import os
from itertools import chain, zip_longest
import pickle


def plyify(r, t, ply_thickness, reverse=False):
    """Fill thickness distribution using plies."""
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
    if reverse:
        return done
    else:
        return [i for i in reversed(done)]


def ply_stack(r, t, t_ply=1.0, reverse=False, subdivisions=5000, material=11):
    """
    Taking in a r,t thickness over radius distribution, split it up in
    plies, note that it rounds down in thickness (if the last ply does not fit
    in the thickness it does not go on), so if you have thick plies,
    make sure that the mm thickness is slightly up from a whole number of plies.
    """
    x = np.linspace(min(r), max(r), subdivisions)
    y = np.interp(x, list(r), t)
    st = plyify(x, y, t_ply, reverse)
    return [[material, t_ply] + i for i in st]


def coreblock(r, t, subdivisions=200, material=11):
    """Make a thickness distribution into a block divisions."""

    x = np.array(sorted(list(np.linspace(min(r), max(r), int(subdivisions))) + list(r)))
    y = np.interp(0.5 * (x[:-1] + x[1:]), list(r), t)

    if type(material) == str and material.startswith(
        "["
    ):  # deal with multi-material cores here, parse the multi-material distribution
        mat = list(zip(*eval(material)))
        mm = np.interp(np.linspace(0, 1, len(x)), mat[1], mat[0])
        mm = [min(mat[0], key=lambda x: abs(i - x)) for i in mm]
    elif type(material) == int or type(eval(material)) == int:
        mm = np.ones(len(x)) * int(material)

    stack = [i for i in list(zip(mm, y, x[:-1], x[1:])) if i[1] != 0]
    return stack


def number_stack(stack, splitstack, key, increment):
    """Create ply numbering for the plies in the stack, depends on the start key and increment as well as the above/below ratio (splitstack)."""
    sp = np.nan_to_num(splitstack.astype(float))
    if sp.sum() != 1.0:
        exit("sum %s not 1." % splitstack)

    ky = np.nan_to_num(key.astype(float)).astype(int)
    inc = np.nan_to_num(increment.astype(float), nan=1).astype(int)

    ss = (len(stack) * sp).round().astype(int)
    st, sb = np.array(range(ss[0])), np.array(range(ss[1]))
    stt = ky[0] + st * inc[0]
    sbt = ky[1] + sb * inc[1]
    seq = np.array(
        [i for i in chain.from_iterable(zip_longest(stt, sbt)) if i is not None]
    ).astype(int)
    return seq


def get_coverage(coverage, parameters, rr):
    # replace the parameter names by array
    for p in parameters:
        if coverage.find(p) != -1:
            coverage = coverage.replace(p, repr(parameters[p]))

    o = []
    for i in coverage:
        out = [
            i,
            np.ones(len(rr)) * coverage[i][0],
            np.ones(len(rr)) * coverage[i][1],
            np.ones(len(rr)) * coverage[i][2],
        ] + [rr]
        o.append(out)

    return o


def lamplan2plies(blade):
    root_radius = blade["planform"]["z"][0][1]
    tip_radius = blade["planform"]["z"][-1][1]

    slabs = blade["laminates"]["slabs"]

    allstacks = []
    parameters = {}  # parameters defined in the laminate plan csv file

    # use a multiple of lamplan length as radius grid to interpolate geometric variables to
    n_s = round(tip_radius * 4)

    rr = np.linspace(root_radius, tip_radius, n_s)

    material_map = {}

    for i in slabs:
        name = i
        material = "glass_biax" if "material" not in slabs[i] else slabs[i]["material"]
        if material not in material_map:
            material_map[material] = len(material_map)
        grid = "lamplan" if "grid" not in slabs[i] else slabs[i]["grid"]
        coverage = (
            slabs[i]["cover"]
            if "cover" in slabs[i]
            else exit("no cover defined for slab %s" % i)
        )
        draping = "plies" if "draping" not in slabs[i] else slabs[i]["draping"]
        splitstack = (
            np.array([1, 0])
            if "splitstack" not in slabs[i]
            else np.array(slabs[i]["splitstack"])
        )
        key = (
            np.array([0, 4000]) if "key" not in slabs[i] else np.array(slabs[i]["key"])
        )
        increment = (
            np.array([1, -1])
            if "increment" not in slabs[i]
            else np.array(slabs[i]["increment"])
        )
        scale = (
            tip_radius - root_radius if "rscale" not in slabs[i] else slabs[i]["rscale"]
        )
        ply_thickness = (
            1.0 if "ply_thickness" not in slabs[i] else slabs[i]["ply_thickness"]
        )
        r, t = np.array(slabs[i]["slab"]).T
        r *= scale
        t *= ply_thickness

        cover = get_coverage(coverage, parameters, rr)
        if draping == "blocks":
            stack = coreblock(r, t, material=material_map[material])
        elif draping == "plies":
            stack = ply_stack(
                r, t, float(ply_thickness), material=material_map[material]
            )

        # assigns keys to the plies in the stack depending on how the stack is split up
        # what increment the ply keys are stacked with (non-1 increment allowing interleaving of plies)
        # print(splitstack, key)
        stack_numbering = number_stack(stack, splitstack, key, increment)
        allstacks.append((name.strip(), grid.strip(), cover, stack_numbering, stack))

    print("material map ", material_map)

    return allstacks


def main():
    parser = argparse.ArgumentParser(
        description="Split up the slab based laminate plan into plies (for laminae) and blocks (for core materials), write to a .pck file for use in draping"
    )
    parser.add_argument("yaml", help="Laminate plan yaml file")
    parser.add_argument("--out", default="__lamplan.pck", help="Output file name")
    args = parser.parse_args()

    blade = yaml.load(open(args.yaml, "r"), Loader=yaml.CLoader)

    allstacks = lamplan2plies(blade)

    of = args.out
    pickle.dump(allstacks, open(of, "wb"))

    print("written plydrape to %s" % of)


if __name__ == "__main__":
    main()
