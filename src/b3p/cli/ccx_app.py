import logging
from functools import partial
from pathlib import Path
import os
import glob
import multiprocessing
import subprocess
from rich.progress import Progress
from rich.logging import RichHandler
from ..ccx import mesh2ccx, ccx2vtu, ccxpost
from ..ccx.failcrit_mesh import compute_failure_for_meshes
from statesman.core.base import Statesman, ManagedFile
from treeparse import cli, command, option

logger = logging.getLogger(__name__)
# logging.basicConfig(handlers=[RichHandler(rich_tracebacks=True)], level=logging.INFO)


class PrepStep(Statesman):
    """Step for preparing CCX input files."""

    dependent_sections = ["mesh"]
    output_files = ["ccx_input.inp"]

    def _execute(self, bondline=False):
        base_prefix = self.workdir / self.config["general"]["prefix"]
        prefix = self.workdir / "fea" / self.config["general"]["prefix"]
        available_meshes = glob.glob(f"{base_prefix}_joined.vtu")
        if bondline:
            bondline_meshes = glob.glob(f"{base_prefix}*_bondline.vtu")
            if bondline_meshes:
                available_meshes = bondline_meshes

        if not available_meshes:
            logger.error("No meshes found, did you build the blade geometry?")
            return

        output_files = mesh2ccx.mesh2ccx(
            available_meshes[-1],
            matmap=str(self.workdir / "material_map.json"),
            out=f"{prefix}_ccx.inp",
            bondline=bondline,
        )
        logger.info(f"Written: {', '.join(output_files)}")


class SolveStep(Statesman):
    """Step for solving CCX problem."""

    input_files = [
        ManagedFile(name="ccx_input.inp", non_empty=True),
    ]
    output_files = ["ccx_output.frd"]
    dependent_sections = ["ccx"]

    def _execute(self, wildcard="", nproc=2, ccxexe="ccx"):
        prefix = self.workdir / "fea" / self.config["general"]["prefix"]
        inpfiles = glob.glob(f"{prefix}*ccx*{wildcard}*.inp")
        inps_to_run = [inp for inp in inpfiles if not check_ccx_run_done(inp)]

        if inps_to_run:
            with multiprocessing.Pool(nproc) as pool:
                with Progress() as progress:
                    task = progress.add_task("Running CCX", total=len(inps_to_run))
                    for inp, success, error_msg in pool.imap_unordered(
                        partial(run_ccx, ccxexe=ccxexe, logger=logger), inps_to_run
                    ):
                        progress.update(task, advance=1)
                        if not success:
                            logger.error(error_msg)


def check_ccx_run_done(inpfile):
    frd_file = inpfile.replace(".inp", ".frd")
    if os.path.exists(frd_file):
        with open(frd_file, "rb") as f:
            f.seek(-5, 2)
            y = f.read()
            if y == b"9999\n":
                logger.info(f"CCX run for {inpfile} is done, found {frd_file}")
                return True
    return False


def run_ccx(inp, ccxexe, logger):
    cmd = [ccxexe, inp.replace(".inp", "")]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return inp, True, None
    except subprocess.CalledProcessError as e:
        error_msg = f"ccx failed for {inp}: {e}"
        return inp, False, error_msg


# CLI code for ccx remains, but integrated with statesman
# ... (rest unchanged)


def run_callback(yml: Path, bondline: bool, buckling: bool):
    state = AppState.get_instance()
    app = CcxApp(state, yml)
    app.ccx(bondline=bondline, buckling=buckling)


def prep_callback(yml: Path, bondline: bool, buckling: bool):
    state = AppState.get_instance()
    app = CcxApp(state, yml)
    app.prep(bondline=bondline, buckling=buckling)


def solve_callback(
    yml: Path, wildcard: str, nproc: int, ccxexe: str, merged_plies: bool
):
    state = AppState.get_instance()
    app = CcxApp(state, yml)
    app.solve(wildcard=wildcard, nproc=nproc, ccxexe=ccxexe, merged_plies=merged_plies)


def post_callback(yml: Path, wildcard: str, nbins: int):
    state = AppState.get_instance()
    app = CcxApp(state, yml)
    app.post(wildcard=wildcard, nbins=nbins)


def plot_callback(yml: Path, plot3d: bool, plot2d: bool):
    state = AppState.get_instance()
    app = CcxApp(state, yml)
    app.plot(plot3d=plot3d, plot2d=plot2d)


def failure_callback(yml: Path):
    state = AppState.get_instance()
    app = CcxApp(state, yml)
    app.failure_criteria()


ccx_cli = cli(
    name="ccx",
    help="Run Calculix operations",
    line_connect=True,
    show_types=True,
    show_defaults=True,
    options=[
        option(
            flags=["--yml", "-y"],
            arg_type=Path,
            required=True,
            help="Path to YAML config file",
        ),
    ],
)

ccx_cli.commands.append(
    command(
        name="run",
        help="Run full Calculix process",
        callback=run_callback,
        arguments=[],
        options=[
            option(
                flags=["--bondline", "-b"],
                is_flag=True,
                default=False,
                arg_type=bool,
                help="Use bondline meshes",
            ),
            option(
                flags=["--buckling", "-k"],
                is_flag=True,
                default=False,
                arg_type=bool,
                help="Enable buckling analysis",
            ),
        ],
    )
)

ccx_cli.commands.append(
    command(
        name="prep",
        help="Prepare CCX input files",
        callback=prep_callback,
        arguments=[],
        options=[
            option(
                flags=["--bondline", "-b"],
                is_flag=True,
                default=False,
                arg_type=bool,
                help="Use bondline meshes",
            ),
            option(
                flags=["--buckling", "-k"],
                is_flag=True,
                default=False,
                arg_type=bool,
                help="Enable buckling analysis",
            ),
        ],
    )
)

ccx_cli.commands.append(
    command(
        name="solve",
        help="Solve CCX problem",
        callback=solve_callback,
        arguments=[],
        options=[
            option(
                flags=["--wildcard", "-w"],
                default="",
                arg_type=str,
                help="Wildcard pattern for input files",
            ),
            option(
                flags=["--nproc", "-p"],
                default=2,
                arg_type=int,
                help="Number of processes",
            ),
            option(
                flags=["--ccxexe", "-c"],
                default="ccx",
                arg_type=str,
                help="Calculix executable",
            ),
            option(
                flags=["--merged-plies", "-m"],
                is_flag=True,
                default=False,
                arg_type=bool,
                help="Only process merged plies",
            ),
        ],
    )
)

ccx_cli.commands.append(
    command(
        name="post",
        help="Postprocess CCX results",
        callback=post_callback,
        arguments=[],
        options=[
            option(
                flags=["--wildcard", "-w"],
                default="",
                arg_type=str,
                help="Wildcard pattern for results",
            ),
            option(
                flags=["--nbins", "-n"],
                default=60,
                arg_type=int,
                help="Number of bins for tabulation",
            ),
        ],
    )
)

ccx_cli.commands.append(
    command(
        name="plot",
        help="Plot CCX results",
        callback=plot_callback,
        arguments=[],
        options=[
            option(
                flags=["--plot3d", "-3"],
                is_flag=True,
                default=True,
                arg_type=bool,
                help="Enable 3D plots",
            ),
            option(
                flags=["--plot2d", "-2"],
                is_flag=True,
                default=True,
                arg_type=bool,
                help="Enable 2D plots",
            ),
        ],
    )
)

ccx_cli.commands.append(
    command(
        name="failure",
        help="Compute failure criteria",
        callback=failure_callback,
        arguments=[],
    )
)
