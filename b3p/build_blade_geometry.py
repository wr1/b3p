#! /usr/bin/env python3

import os
from copy import deepcopy as dc
import argparse
import b3p.blade
import b3p.loft_utils
import yaml
import numpy as np


def run_loft(configfile, verbose=False):
    """
    Perform blade 3d model generation based on yaml input file.
    """
    config = yaml.load(open(configfile, "r"), Loader=yaml.CLoader)

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
        barrel_length=2,
        interpolate_method=0,
        flatten_lw=False,
        offset_optimal=True,
        offset_clamp_points=[-0.001, 0.04, 0.2, 0.45, 0.7, 0.8, 0.9, 0.97],
    )

    wdp = os.path.join(config["general"]["workdir"], config["general"]["prefix"])

    if not os.path.isdir(config["general"]["workdir"]):
        os.makedirs(config["general"]["workdir"])

    blade.mesh("%s.stl" % wdp)
    blade.dump("%s.pck" % (wdp), z_rotation=0)
    blade.export_variables("%s.var" % (wdp))
    blade.export_xfoil(
        prefix=os.path.join(
            config["general"]["workdir"],
            "airfoil_out",
            "xs_%s" % config["general"]["prefix"],
        )
    )
    blade.plot("%s" % wdp, fname="%s.png" % wdp)
    n_sections = 50
    blade.to_table("%s_sca_%i" % (wdp, n_sections), np.linspace(0, 1, n_sections))


def main():
    parser = argparse.ArgumentParser(
        description="Generate a blade from a configuration file"
    )
    parser.add_argument("yaml", help="yaml blade definition file")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    run_loft(args.yaml, args.verbose)


if __name__ == "__main__":
    main()
