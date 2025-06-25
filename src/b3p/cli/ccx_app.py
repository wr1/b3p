import logging
from functools import partial
from pathlib import Path
import os
import glob
import multiprocessing
import subprocess
from rich.progress import Progress  # Replace tqdm with rich progress
from rich.logging import RichHandler  # Add rich log formatting
from b3p.ccx import mesh2ccx, ccx2vtu, ccxpost
from b3p.ccx.failcrit_mesh import compute_failure_for_meshes
from b3p.cli.app_state import AppState

logger = logging.getLogger(__name__)
logging.basicConfig(handlers=[RichHandler(rich_tracebacks=True)], level=logging.INFO)
# Configure rich for logging


def run_ccx(inp, ccxexe, logger):
    """Run CalculiX (ccx) on a given input file, capturing stdout/stderr."""
    cmd = [ccxexe, inp.replace(".inp", "")]
    logger.debug(f"Running command: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.debug(f"CCX output for {inp}:\n{result.stdout}")
        return inp, True, None
    except subprocess.CalledProcessError as e:
        error_msg = (
            f"ccx failed for {inp}: {e}\nStdout:\n{e.stdout}\nStderr:\n{e.stderr}"
        )
        logger.error(error_msg)
        return inp, False, error_msg


def check_ccx_run_done(inpfile):
    """Check if the ccx run is done by looking for the .frd file."""
    frd_file = inpfile.replace(".inp", ".frd")
    if os.path.exists(frd_file):
        with open(frd_file, "rb") as f:
            f.seek(-5, 2)
            y = f.read()
            if y == b"9999\n":
                return True
    return False


class CcxApp:
    def __init__(self, state: AppState, yml: Path):
        self.state = state
        self.yml = yml
        self.dir = "fea"
        self.state.load_yaml(self.yml)

    def ccx(self, **kwargs):
        self.prep(**kwargs)
        self.solve(**kwargs)
        self.post(**kwargs)
        self.plot(**kwargs)

    def prep(self, bondline=False, **kwargs):
        base_prefix = self.state.get_prefix("drape")
        prefix = self.state.get_prefix(self.dir)
        available_meshes = glob.glob(f"{base_prefix}_joined.vtu")
        logger.info(f" available meshes {available_meshes}")
        if bondline:
            bondline_meshes = glob.glob(f"{base_prefix}*_bondline.vtu")
            if bondline_meshes:
                available_meshes = bondline_meshes

        logger.info(f"Available meshes: {available_meshes}")
        if not available_meshes:
            logger.error("No meshes found, did you build the blade geometry?")
            return

        output_files = mesh2ccx.mesh2ccx(
            available_meshes[-1],
            matmap=str(Path(base_prefix).parent / "material_map.json"),
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
        bondline=False,
        **kwargs,
    ):
        self.state.load_yaml(self.yml)
        prefix = self.state.get_prefix(self.dir)
        if inpfiles is None:
            inpfiles = glob.glob(f"{prefix}*ccx*{wildcard}*.inp")
        inps = [inp for inp in inpfiles]

        inps_to_run = [inp for inp in inps if not check_ccx_run_done(inp)]

        if not inps:
            logger.error(f"No input files found matching {prefix}*{wildcard}*inp")
            return
        else:
            logger.info(
                f"Found {len(inps)} input files matching {prefix}*{wildcard}*inp"
            )

        if inps_to_run:
            with multiprocessing.Pool(nproc) as pool:
                with Progress() as progress:  # Replace tqdm with rich Progress
                    task = progress.add_task("Running CCX", total=len(inps_to_run))
                    for inp, success, error_msg in pool.imap_unordered(
                        partial(run_ccx, ccxexe=ccxexe, logger=logger), inps_to_run
                    ):
                        progress.update(task, advance=1)
                        if not success:
                            logger.error(error_msg)

    def post(self, wildcard="", nbins=60, bondline=False, **kwargs):
        self.state.load_yaml(self.yml)
        prefix = self.state.get_prefix(self.dir)
        ccxpost = ccx2vtu.ccx2vtu(prefix, wildcard=wildcard)
        ccxpost.load_grids()
        ccxpost.tabulate(nbins)

    def failure_criteria(self, **kwargs):
        self.state.load_yaml(self.yml)
        prefix = self.state.get_prefix(self.dir)
        puck_config = self.state.dct.get("damage", {})
        vtus = [i for i in prefix.parent.glob("*ccx*vtu") if str(i).find("fail") == -1]
        logger.info(f"VTU files found for failure criteria: {vtus}")

        compute_failure_for_meshes(vtus, puck_config)

    def plot(self, plot3d=True, plot2d=True, bondline=False, **kwargs):
        self.state.load_yaml(self.yml)
        plotter = ccxpost.plot_ccx(self.state.get_workdir())

        print(f"plotting 2d {plot2d} and 3d {plot3d}")
        if plot3d:
            plotter.plot3d()
        if plot2d:
            plotter.plot2d()
