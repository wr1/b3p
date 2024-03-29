#! /usr/bin/env python

from ruamel import yaml
import numpy as np
from ruamel.yaml import YAML
import re
import fire
import os


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
    """Load an airfoil from an xfoil formated text file.
    Check if the first line are coordinates or airfoil name."""
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


def load_airfoils(afs, prefix=""):
    """Load a list of airfoils from a dictionary of xfoil formated text files."""
    dct = {}
    for i in afs:
        if "xy" in afs[i]:
            return afs
        name, xy = load_airfoil(os.path.join(prefix, afs[i]))
        dct[i] = {"xy": flowlist(xy)}
        dct[i]["name"] = name
        dct[i]["path"] = afs[i]
        print(f"\t** imported airfoil {name} at thickness {i}")

    return dct


def yaml_make_portable(yaml_file, safe=False):
    """Import a yaml file and all linked files, and write to a _portable version of itself.

    :param yaml_file: path to yaml file
    :param safe: if True, use safe loader and basic data types (dict, list, str, int, float)
    :return: None"""
    print(f"** Loading yaml file {yaml_file} and loading linked files")
    prefix = os.path.dirname(yaml_file)

    yaml = YAML(typ="safe") if safe else YAML(typ="rt")
    d = yaml.load(open(yaml_file, "r"))

    d["aero"]["airfoils"] = load_airfoils(d["aero"]["airfoils"], prefix=prefix)

    if type(d["materials"]) == str:
        print(f'\t** loading materials from: {d["materials"]}')
        d["materials"] = yaml.load(open(os.path.join(prefix, d["materials"]), "r"))

    if d["general"]["workdir"].find("portable") == -1:
        d["general"]["workdir"] += "_portable"

    return d


def save_yaml_portable(yaml_file):
    d = yaml_make_portable(yaml_file)

    y = YAML()
    of = yaml_file.replace(".yml", "_portable.yml").replace(".yaml", "_portable.yml")
    with open(of, "w") as f:
        y.dump(d, f)
    print(f"\twritten to: {of}")


def main():
    fire.Fire(save_yaml_portable)


if __name__ == "__main__":
    main()
