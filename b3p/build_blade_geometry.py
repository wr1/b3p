#! /usr/bin/env python3

import os
from copy import deepcopy as dc
import fire
import b3p.blade
import b3p.loft_utils
from ruamel import yaml
import numpy as np


def build_blade_geometry(config, verbose=False, xfoil=False):
    """
    Perform blade 3d model generation based on b3p dictionary.

    :param config: b3p dictionary
    """
    pln = config["planform"]

    blade = b3p.blade.blade(
        pln["chord"],
        pln["thickness"],
        pln["twist"],
        dc(pln["dx"]),
        pln["dy"],
        pln["z"],
        config["aero"]["airfoils"],
        chordwise_sampling=b3p.loft_utils.optspace(config["planform"]["npchord"]),
        np_spanwise=config["planform"]["npspan"],
    )

    wdp = os.path.join(config["general"]["workdir"], config["general"]["prefix"])

    if not os.path.isdir(config["general"]["workdir"]):
        os.makedirs(config["general"]["workdir"])

    blade.mesh(f"{wdp}.vtp")
    blade.dump(f"{wdp}.pck", z_rotation=0)
    blade.export_variables(f"{wdp}.var")
    if xfoil:
        blade.export_xfoil(
            prefix=os.path.join(
                config["general"]["workdir"],
                "airfoil_out",
                f'xs_{config["general"]["prefix"]}',
            )
        )
    blade.plot(f"{wdp}", fname=f"{wdp}.png")
    n_sections = 50
    blade.to_table(np.linspace(0, 1, n_sections), "%s_sca_%i" % (wdp, n_sections))
    return blade


def run_blade_geometry(yamlfile, verbose=False):
    """Build a blade geometry from a yaml file"""
    config = yaml.load(open(yamlfile, "r"), Loader=yaml.CLoader)
    _ = build_blade_geometry(config, verbose=verbose)


def main():
    fire.Fire(run_blade_geometry)
