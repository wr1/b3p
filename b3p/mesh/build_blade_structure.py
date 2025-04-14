#! /usr/bin/env python3

from b3p.mesh import mesh_from_loft
from b3p.mesh import webs
import numpy as np


def build_blade_structure(config, prefix=None):
    pckfile = f"{prefix}.pck"
    base_vtp = f"{prefix}_base.vtp"
    shell_vtp = f"{prefix}_shell.vtp"

    radii = eval(str(config["mesh"]["radii"]))
    config_webs = config["mesh"]["webs"]
    panel_mesh_scale = (
        config["mesh"]["panel_mesh_scale"]
        if "panel_mesh_scale" in config["mesh"]
        else []
    )

    # build the plain blade mesh
    mesh_from_loft.build_mesh(pckfile, radii, [], [], prefix, outfile=base_vtp)

    # compute the locations where the webs intersect with the blade mesh
    web_shell_intersections = webs.build_webs(base_vtp, config_webs, prefix=prefix)

    # rejoggle the datums format from the yaml file into dict
    if "coordinates" in config["mesh"]:
        ad = config["mesh"]["coordinates"]
        added_datums = dict(
            [(i, [ad[i]["base"]] + list(zip(*ad[i]["points"]))) for i in ad]
        )
    else:
        added_datums = {}

    # build the blade mesh with webs
    return mesh_from_loft.build_mesh(
        pckfile,
        radii,
        config_webs,
        web_shell_intersections,
        prefix,
        config["mesh"]["n_web_points"],
        config["mesh"]["n_chordwise_points"],
        outfile=shell_vtp,
        added_datums=added_datums,
        panel_mesh_scale=panel_mesh_scale,
    )
