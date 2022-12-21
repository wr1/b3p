#! /usr/bin/env python3

import os
from copy import deepcopy as dc
from b3p import mesh_from_loft
from b3p import webs
import yaml
import argparse
import numpy as np


def build_rectangle_blade_mesh_with_webs(configfile):
    config = yaml.load(open(configfile, "r"), Loader=yaml.CLoader)
    wdp = os.path.join(config["general"]["workdir"], config["general"]["prefix"])
    pckfile = wdp + ".pck"
    base_vtp = "%s_base.vtp" % wdp
    web_vtp = "%s_web.vtp" % wdp
    radii = eval(str(config["mesh"]["radii"]))
    config_webs = config["mesh"]["webs"]
    panel_mesh_scale = (
        config["mesh"]["panel_mesh_scale"]
        if "panel_mesh_scale" in config["mesh"]
        else []
    )

    # build the plain blade mesh
    mesh_from_loft.build_mesh(pckfile, radii, [], [], wdp, outfile=base_vtp)

    # compute the locations where the webs intersect with the blade mesh
    web_shell_intersections = webs.build_webs(base_vtp, config_webs, prefix=wdp)

    # rejoggle the datums format from the yaml file into dict
    ad = config["mesh"]["datums"]
    added_datums = dict(
        [(i, [ad[i]["base"]] + list(zip(*ad[i]["points"]))) for i in ad]
    )

    # remesh the blade, but now making sure there are nodes at the web intersections
    mesh_from_loft.build_mesh(
        pckfile,
        radii,
        config_webs,
        web_shell_intersections,
        wdp,
        config["mesh"]["n_web_points"],
        config["mesh"]["n_chordwise_points"],
        outfile=web_vtp,
        added_datums=added_datums,
        panel_mesh_scale=panel_mesh_scale,
    )


def main():
    parser = argparse.ArgumentParser(
        "Build a rectangle blade mesh with shared nodes where the webs link up with the shell"
    )
    parser.add_argument("yaml", help="config file ")
    args = parser.parse_args()
    build_rectangle_blade_mesh_with_webs(args.yaml)


if __name__ == "__main__":
    main()
