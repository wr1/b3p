#! /usr/bin/env python3

import os
from copy import deepcopy as dc
import b3p.geometry.blade
import b3p.geometry.loft_utils
import numpy as np


def build_blade_geometry(config, prefix, xfoil=True):
    """Perform blade 3d model generation based on b3p dictionary."""
    pln = config["planform"]
    blade = b3p.geometry.blade.blade(
        pln["chord"],
        pln["thickness"],
        pln["twist"],
        dc(pln["dx"]),
        pln["dy"],
        pln["z"],
        config["aero"]["airfoils"],
        chordwise_sampling=b3p.geometry.loft_utils.optspace(
            config["planform"]["npchord"]
        ),
        np_spanwise=config["planform"]["npspan"],
    )

    blade.mesh(f"{prefix}.vtp")
    blade.dump(f"{prefix}.pck", z_rotation=0)
    blade.export_variables(f"{prefix}_variables.json")
    if xfoil:
        blade.export_xfoil(
            prefix=os.path.join(
                prefix.parent,
                "airfoil_out",
                f"xs_{config['general']['prefix']}",
            )
        )
    blade.plot(f"{prefix}", fname=f"{prefix}.png")
    n_sections = 50
    blade.to_table(np.linspace(0, 1, n_sections), "%s_sca_%i" % (prefix, n_sections))

    return blade
