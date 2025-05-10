import logging
from functools import partial
from pathlib import Path
import os
import glob
import multiprocessing
import subprocess
from b3p.ccx import mesh2ccx, ccx2vtu, ccxpost

logger = logging.getLogger(__name__)


def run_ccx(inp, ccxexe, logger):
    """
    Run CalculiX (ccx) on a given input file.

    Args:
        inp (str): Path to the input file.
        ccxexe (str): Path or name of the ccx executable.
        logger (logging.Logger): Logger instance for logging messages.
    """
    cmd = [ccxexe, inp.replace(".inp", "")]
    logger.info(f"Running command: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"ccx failed for {inp}: {e}")


class CcxApp:
    def __init__(self, state, yml: Path):
        self.state = state
        self.yml = yml
        self.dir = "fea"

    def ccx(self, **kwargs):
        self.prep(**kwargs)
        self.solve(**kwargs)
        self.post(**kwargs)
        self.plot(**kwargs)

    def prep(self, bondline=False, **kwargs):
        dct = self.state.load_yaml(self.yml)
        base_prefix = self.state.get_prefix("drape")
        prefix = self.state.get_prefix(self.dir)
        available_meshes = glob.glob(f"{base_prefix}_joined.vtu")
        logger.debug(f"Available meshes: {prefix}")
        if bondline:
            bondline_meshes = glob.glob(f"{prefix}*_bondline.vtu")
            if bondline_meshes:
                available_meshes = bondline_meshes

        logger.info(f"Available meshes: {available_meshes}")
        if not available_meshes:
            logger.error("No meshes found, did you build the blade geometry?")
            return

        output_files = mesh2ccx.mesh2ccx(
            available_meshes[-1],
            matmap=os.path.join(os.path.dirname(base_prefix), "material_map.json"),
            out=f"{prefix}_ccx.inp",
            bondline=bondline,
            **{k: v for k, v in kwargs.items() if k != "bondline"},
        )
        logger.info(f"Written: {', '.join(output_files)}")

    def solve(
        self,
        wildcard="",
        nproc=2,
        ccxexe="ccx",
        inpfiles=None,
        merged_plies=False,
        bondline=False,
        **kwargs,
    ):
        dct = self.state.load_yaml(self.yml)
        prefix = self.state.get_prefix(self.dir)
        if inpfiles is None:
            inpfiles = glob.glob(f"{prefix}*ccx*{wildcard}*.inp")
        inps = [inp for inp in inpfiles]
        if merged_plies:
            inps = [inp for inp in inps if "_mp_" in inp]
        if not inps:
            logger.error(f"No input files found matching {prefix}*{wildcard}*inp")
            return
        inps_to_run = []
        for inp in inps:
            frd_file = inp.replace(".inp", ".frd")
            if not os.path.exists(frd_file):
                inps_to_run.append(inp)
            else:
                with open(frd_file, "rb") as f:
                    f.seek(-5, 2)
                    y = f.read()
                    if y != b"9999\n":
                        inps_to_run.append(inp)
                    else:
                        logger.info(f"Skipping {inp} because {frd_file} exists")
        logger.info(f"Running ccx on: {inps_to_run} using {wildcard}")

        with multiprocessing.Pool(nproc) as p:
            p.map(partial(run_ccx, ccxexe=ccxexe, logger=logger), inps_to_run)

    def post(self, wildcard="", nbins=60, bondline=False, **kwargs):
        dct = self.state.load_yaml(self.yml)
        prefix = self.state.get_prefix(self.dir)
        ccxpost = ccx2vtu.ccx2vtu(prefix, wildcard=wildcard)
        ccxpost.load_grids()
        ccxpost.tabulate(nbins)

    def plot(self, wildcard="", plot3d=True, plot2d=True, bondline=False, **kwargs):
        dct = self.state.load_yaml(self.yml)
        plotter = ccxpost.plot_ccx(dct["general"]["workdir"], wildcard=wildcard)
        if plot3d:
            plotter.plot3d()
        if plot2d:
            plotter.plot2d(wildcard=wildcard)
