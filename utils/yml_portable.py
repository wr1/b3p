#! /usr/bin/env python

import numpy as np
import re
import fire
import os

# import yaml
from ruamel.yaml import YAML


# def deep_convert_to_float(data):
#     """
#     Recursively convert all string representations of numbers in a nested data structure
#     (including dictionaries and lists) to floats.
#     """
#     if isinstance(data, dict):
#         return {k: deep_convert_to_float(v) for k, v in data.items()}
#     elif isinstance(data, list):
#         return [deep_convert_to_float(item) for item in data]
#     elif isinstance(data, str):
#         try:
#             return float(data)
#         except ValueError:
#             return data  # Return the original string if conversion is not possible
#     else:
#         return data  # Return the data as is if it's not a list, dict, or string


# # Define a custom representer for lists to make them inline
# def custom_list_representer(dumper, data):
#     return dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)


# # Register the custom list representer with the PyYAML dumper
# yaml.add_representer(list, custom_list_representer)


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
        print(f"\t** imported airfoil {name} at thickness {i}")

    return dct


def yaml_make_portable(yaml_file, safe=False):
    """Import a yaml file and all linked files, and write to a _portable version of itself.

    :param yaml_file: path to yaml file
    :param safe: if True, use safe loader and basic data types (dict, list, str, int, float)
    :return: None"""
    print(f"** Loading yaml file {yaml_file} and loading linked files")
    prefix = os.path.dirname(yaml_file)

    yaml = YAML()  # typ="safe")
    # typ="safe")  # YAML(typ="safe") if safe else YAML(typ="rt")
    d = yaml.load(open(yaml_file, "r"))
    # yaml.load(open(yaml_file, "r"), Loader=yaml.CLoader)

    d["aero"]["airfoils"] = load_airfoils(d["aero"]["airfoils"], prefix=prefix)

    if type(d["materials"]) == str:
        print(f'\t** loading materials from: {d["materials"]}')
        d["materials"] = yaml.load(open(os.path.join(prefix, d["materials"]), "r"))

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
        print(f"\twritten to: {of}")


def save_yaml_portable(yaml_file):
    d = yaml_make_portable(yaml_file)
    of = yaml_file.replace(".yml", "_portable.yml").replace(".yaml", "_portable.yml")
    save_yaml(of, d)


def main():
    fire.Fire(save_yaml_portable)


if __name__ == "__main__":
    main()
