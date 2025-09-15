# Runner class for ccblade_app.

from pathlib import Path
import numpy as np
import pandas as pd
from ccblade.ccblade import CCBlade
import os
import logging
from typing import List, Tuple, Dict, Optional, Any

from .utils import load_polar, interpolate_polars
from .optimizer import controloptimize

logger = logging.getLogger(__name__)


class ccblade_run:
    """Run CCBlade analysis on a blade."""

    def __init__(self, config: dict, yml_dir: Path):
        """Initialize with a config dict and YAML directory."""
        self.config = config
        self.yml_dir = yml_dir

        # Resolve workdir relative to YAML directory
        workdir_path = Path(self.config["general"]["workdir"])
        if not workdir_path.is_absolute():
            workdir_path = self.yml_dir / workdir_path
        self.workdir = workdir_path.resolve() / "mesh"
        self.workdir.mkdir(parents=True, exist_ok=True)  # Ensure mesh directory exists

        bem = self.config["aero"]["bem"]
        self.prefix = self.workdir / self.config["general"]["prefix"]

        if bem["polars"] is None:
            exit("no polars in blade file")
        if not os.path.isfile(self.prefix.with_suffix(".pck")):
            exit("blade not built yet, run b3p build <blade.yml> first")
        plf = pd.read_csv(
            self.prefix.with_name(self.prefix.stem + "_sca_50.csv"), sep=";"
        )
        plrs = sorted(
            [(i['key'], load_polar(yml_dir / Path(i['file']))) for i in bem["polars"]],
            reverse=True,
        )
        iplr = interpolate_polars(
            plrs, plf.relative_thickness, of=self.workdir.parent / "polars.png"
        )
        rhub, rtip = plf.z.iloc[0], plf.z.iloc[-1]
        self.rotor = CCBlade(
            plf.z - plf.z[0],
            plf.chord,
            plf.twist,
            iplr,
            rhub,
            rtip,
            B=bem["B"],
            rho=bem["rho"],
            mu=bem["mu"],
            precone=bem["precone"],
            tilt=bem["tilt"],
            yaw=bem["yaw"],
            shearExp=bem["shearExp"],
            hubHt=bem["hubHt"],
            derivatives=False,
        )
        logger.info(f"Rotor from {rhub} to {rtip}")
        self.copt = controloptimize(
            self.rotor,
            bem["max_tipspeed"],
            rtip,
            bem["rated_power"],
            uinf=np.array(bem["uinf"]),
            workdir=self.workdir,
        )

    def run(self) -> None:
        """Execute the CCBlade analysis."""
        self.copt.control_opt_below_rated()
        output = self.copt.control_opt_above_rated()
        self.copt.compute_bladeloads()
        del output["W"]
        df = pd.DataFrame(output).dropna()
        df.to_csv(self.workdir.parent / "ccblade_output.csv", sep=";")
