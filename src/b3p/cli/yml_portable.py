#! /usr/bin/env python

import numpy as np
import re
import os
import logging
from ruamel.yaml import YAML

logger = logging.getLogger(__name__)


def load_airfoil(af):
    """Load an airfoil from an xfoil formatted text file."""
    if not os.path.exists(af):
        logger.warning(f"Airfoil file {af} not found, keeping as reference")
        return None, None
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
    """Load a list of airfoils from a dictionary of xfoil formatted text files."""
    dct = {}
    for i in afs:
        if "xy" in afs[i]:
            dct[i] = afs[i]
            continue
        af_path = os.path.join(prefix, afs[i])
        name, xy = load_airfoil(af_path)
        if name is None:
            dct[i] = {"path": afs[i]}
        else:
            dct[i] = {"xy": xy, "name": name, "path": afs[i]}
            logger.info(f"Imported airfoil {name} at thickness {i}")
    return dct


def yaml_make_portable(yaml_file):
    """Import a yaml file and linked files, skipping missing files."""
    logger.info(f"Loading yaml file {yaml_file}")
    prefix = os.path.dirname(yaml_file)

    yaml = YAML()
    d = yaml.load(open(yaml_file, "r"))

    if "aero" in d and "airfoils" in d["aero"]:
        d["aero"]["airfoils"] = load_airfoils(d["aero"]["airfoils"], prefix=prefix)

    subsections = ["materials", "loads", "laminates"]
    for s in subsections:
        if s in d and isinstance(d[s], str) and d[s].endswith(".yml"):
            sub_path = os.path.join(prefix, d[s])
            if os.path.exists(sub_path):
                logger.info(f"Loading {s} from: {sub_path}")
                d[s] = yaml.load(open(sub_path, "r"))
            else:
                logger.warning(f"File {sub_path} not found, keeping as reference")
                d[s] = {"path": d[s]}

    return d


def save_yaml(of, dct):
    """Save dictionary to yaml file."""
    yaml = YAML()
    yaml.default_flow_style = True
    with open(of, "w") as f:
        yaml.dump(dct, f)
        logger.info(f"Written to: {of}")


def save_yaml_portable(yaml_file):
    """Save portable yaml file."""
    d = yaml_make_portable(yaml_file)
    of = os.path.join(os.path.dirname(yaml_file), "portable.yml")
    save_yaml(of, d)
