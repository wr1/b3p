#! /usr/bin/env python

import numpy as np
import re
import os
import logging
from ruamel.yaml import YAML

logger = logging.getLogger(__name__)


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
        dct[i] = {"xy": xy}
        dct[i]["name"] = name
        dct[i]["path"] = afs[i]
        logger.info(f"imported airfoil {name} at thickness {i}")

    return dct


def yaml_make_portable(yaml_file):
    """Import a yaml file and all linked files, and write to a _portable version of itself.

    :param yaml_file: path to yaml file
    :param safe: if True, use safe loader and basic data types (dict, list, str, int, float)
    :return: None"""
    logger.info(f"Loading yaml file {yaml_file} and loading linked files")
    prefix = os.path.dirname(yaml_file)

    yaml = YAML()
    # typ="safe")
    # typ="safe")  # YAML(typ="safe") if safe else YAML(typ="rt")
    d = yaml.load(open(yaml_file, "r"))
    # yaml.load(open(yaml_file, "r"), Loader=yaml.CLoader)

    d["aero"]["airfoils"] = load_airfoils(d["aero"]["airfoils"], prefix=prefix)

    subsections = ["materials", "loads", "laminates"]

    for s in subsections:
        if type(d[s]) == str and d[s].find(".yml") != -1:
            logger.info(f"loading {s} from: {d[s]}")
            d[s] = yaml.load(open(os.path.join(prefix, d[s]), "r"))

    if d["general"]["workdir"].find("portable") == -1:
        d["general"]["workdir"] += "_portable"
    return d  # deep_convert_to_float(d)


# Alternatively, define a custom representation for lists to ensure all are inline
def represent_list(dumper, data):
    return dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)


def save_yaml(of, dct):
    yaml = YAML()
    yaml.default_flow_style = True
    yaml.representer.add_representer(list, represent_list)
    with open(of, "w") as f:
        yaml.dump(dct, f)
        logger.info(f"written to: {of}")


def save_yaml_portable(yaml_file):
    d = yaml_make_portable(yaml_file)
    of = yaml_file.replace(".yml", "_portable.yml").replace(".yaml", "_portable.yml")
    save_yaml(of, d)
