#! /usr/bin/env python

# from ruamel import YAML as YAML
from ruamel.yaml import YAML
import argparse
from matplotlib import pyplot as plt
import numpy as np


def plot_planform(bl, prefix):
    odct = {}
    absolute_thickness = np.array(bl["chord"]["values"]) * np.array(
        bl["rthick"]["values"]
    )
    fig, ax = plt.subplots(5, 2, figsize=(10, 12))
    c = 0
    bl["twist"]["values"] = np.degrees(bl["twist"]["values"]).tolist()
    for i in bl:
        if "grid" in bl[i] and "values" in bl[i]:
            ax[c][0].plot(bl[i]["grid"], bl[i]["values"])
            ax[c][0].set_title(i)
            odct[i] = [list(j) for j in zip(bl[i]["grid"], bl[i]["values"])]
            c += 1

    ax[c][0].plot(bl["chord"]["grid"], absolute_thickness)

    for c, i in enumerate(bl["reference_axis"]):
        ax[c][1].plot(
            bl["reference_axis"][i]["grid"], bl["reference_axis"][i]["values"]
        )
        ax[c][1].set_title(i)
        odct[("d" if i in ["x", "y"] else "") + i] = [
            list(j)
            for j in zip(
                bl["reference_axis"][i]["grid"], bl["reference_axis"][i]["values"]
            )
        ]
    of = f"{prefix}_planform"
    fig.savefig(of)
    print(f"written plot to {of}")

    odct["thickness"] = odct["rthick"]
    del odct["rthick"]
    return odct


def export_airfoils(af, prefix):
    x, ax = plt.subplots(1, 1, figsize=(12, 12))

    af_dict = {}
    for i in af:
        fn = f'airfoils/af_{i["name"]}.dat'
        buf = "%s\n" % i["name"]
        for j in zip(i["coordinates"]["x"], i["coordinates"]["y"]):
            buf += "%s %s\n" % j
        open(fn, "w").write(buf)

        ax.plot(i["coordinates"]["x"], i["coordinates"]["y"], label=fn)

        af_dict[i["relative_thickness"]] = fn
    ax.legend(loc="best")
    plt.savefig(f"{prefix}_airfoils.png")

    return af_dict


def add_rthick_if_absent(x):
    shp = x["components"]["blade"]["outer_shape_bem"]

    if "rthick" not in shp:
        af = x["airfoils"]

        thick = dict([(i["name"], i["relative_thickness"]) for i in af])

        rthick = [thick[i] for i in shp["airfoil_position"]["labels"]]

        shp["rthick"] = {
            "grid": shp["chord"]["grid"],
            "values": np.interp(
                shp["chord"]["grid"], shp["airfoil_position"]["grid"], rthick
            ).tolist(),
        }

    return x


def plot_laminates(model, prefix):
    lam = model["components"]["blade"]["internal_structure_2d_fem"]

    fig, ax = plt.subplots(1, 1, figsize=(12, 12))

    for i in lam["layers"]:
        ax.plot(i["thickness"]["grid"], i["thickness"]["values"], label=i["name"])

    ax.legend(loc="best")

    fig.savefig(f"{prefix}_lams")


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="read a windio file and convert to b3p input, plot airfoils, planform and composites"
    )
    p.add_argument("yaml", help="WindIO yaml input file")
    p.add_argument("template", help="b3p template yaml file")
    p.add_argument("-o", default="out")
    args = p.parse_args()

    yaml = YAML(typ="rt")
    template = yaml.load(open(args.template, "r"))
    # , Loader=yaml.RoundTripLoader)

    x = yaml.load(open(args.yaml, "r"))

    x = add_rthick_if_absent(x)

    out_var = plot_planform(x["components"]["blade"]["outer_shape_bem"], args.o)

    af_dict = export_airfoils(x["airfoils"], args.o)

    for i in out_var:
        template["planform"][i] = out_var[i]
        # print(i, out_var[i])

    template["aero"]["airfoils"] = af_dict

    #
    plot_laminates(x, args.o)

    yaml.dump(
        dict(template),
        open(f"{args.o}.yml", "w"),
        # Dumper=yaml.RoundTripDumper,
        # default_flow_style=None,
    )
