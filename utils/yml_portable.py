#! /usr/bin/env python

from ruamel import yaml
import argparse
import numpy as np
from ruamel.yaml import YAML
import re


def flowlist(listoflists):
    """Set a list of lists to flow style in ruaeml.yaml.

    :param listoflists: list of lists of floats
    :return: ruamel.yaml sequence
    """
    y = YAML()
    olist = []
    for i in listoflists:
        seq = y.seq([round(j, 4) for j in i])
        seq.fa.set_flow_style()
        olist.append(seq)
    oseq = y.seq(olist)
    oseq.fa.set_flow_style()
    return oseq


def load_airfoil(af):
    with open(af, "r") as f:
        name = f.readline().strip()
        floats = re.findall(r"\d+\.\d+", name)
        if len(floats) == 2:
            offset = 0
            name = af.split("/")[-1].split(".")[0]
        elif name == "":
            offset = 1
            name = af.split("/")[-1].split(".")[0]
        else:
            offset = 1

    return name, np.loadtxt(af, skiprows=offset).tolist()


def load_airfoils(afs):
    dct = {}
    for i in afs:
        name, xy = load_airfoil(afs[i])
        dct[i] = {"xy": flowlist(xy)}
        dct[i]["name"] = name
        dct[i]["path"] = afs[i]

    return dct


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("yaml_file", help="yaml file to load")
    args = parser.parse_args()

    y = YAML()

    d = yaml.round_trip_load(open(args.yaml_file, "r"), preserve_quotes=True)

    d["aero"]["airfoils"] = load_airfoils(d["aero"]["airfoils"])

    d["materials"] = yaml.round_trip_load(
        open(d["materials"], "r"), preserve_quotes=True
    )

    of = args.yaml_file.replace(".yml", "_portable.yml")
    with open(of, "w") as f:
        y.dump(d, f)

    print(f"written to: {of}")
