import argparse
from pathlib import Path
from b3p.cli.app_state import AppState
from b3p.cli.build_app import BuildApp
from b3p.cli.ccx_app import CcxApp
from b3p.cli.two_d_app import TwoDApp
from b3p.cli.ccblade_app import CCBladeApp
from b3p.cli.clean_app import CleanApp
import logging
import sys

from rich.logging import RichHandler  # Add rich log formatting


logging.basicConfig(handlers=[RichHandler(rich_tracebacks=True)], level=logging.INFO)

logger = logging.getLogger(__name__)
# File handler for output.log
file_handler = logging.FileHandler("output.log")
file_handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
logger.addHandler(file_handler)


def main():
    parser = argparse.ArgumentParser(description="Blade Design CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Build command (top-level)
    build_parser = subparsers.add_parser("build", help="Build the full blade model")
    build_parser.add_argument("yml", type=Path, help="Path to YAML config file")
    build_parser.add_argument(
        "-n",
        "--no-bondline",
        action="store_false",
        dest="bondline",
        help="Exclude bondline",
    )
    build_subparsers = build_parser.add_subparsers(dest="subcommand", required=False)

    build_geometry_parser = build_subparsers.add_parser(
        "geometry", help="Build blade geometry"
    )
    build_geometry_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )

    build_mesh_parser = build_subparsers.add_parser("mesh", help="Mesh blade structure")
    build_mesh_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )

    build_drape_parser = build_subparsers.add_parser(
        "drape", help="Drape plies onto mesh"
    )
    build_drape_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )
    build_drape_parser.add_argument(
        "-b", "--bondline", action="store_true", help="Add bondline to mesh"
    )

    build_mass_parser = build_subparsers.add_parser("mass", help="Calculate blade mass")
    build_mass_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )

    build_apply_loads_parser = build_subparsers.add_parser(
        "apply-loads", help="Apply loads to mesh"
    )
    build_apply_loads_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )

    # CCX subcommands
    ccx_parser = subparsers.add_parser("ccx", help="Run Calculix operations")
    ccx_parser.add_argument("yml", type=Path, help="Path to YAML config file")
    ccx_parser.add_argument(
        "-b", "--bondline", action="store_true", help="Use bondline meshes"
    )  # Note: Using -b, assuming no conflict; adjust if needed
    ccx_subparsers = ccx_parser.add_subparsers(dest="subcommand", required=False)

    ccx_ccx_parser = ccx_subparsers.add_parser("ccx", help="Run full Calculix process")
    ccx_ccx_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )
    ccx_ccx_parser.add_argument(
        "-B",
        "--bondline_sub",
        action="store_true",
        dest="bondline",
        help="Use bondline meshes",
    )

    ccx_prep_parser = ccx_subparsers.add_parser("prep", help="Prepare CCX input files")
    ccx_prep_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )
    ccx_prep_parser.add_argument(
        "-B",
        "--bondline_sub",
        action="store_true",
        dest="bondline",
        help="Use bondline meshes",
    )

    ccx_solve_parser = ccx_subparsers.add_parser("solve", help="Solve CCX problem")
    ccx_solve_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )
    ccx_solve_parser.add_argument(
        "-w", "--wildcard", default="", help="Wildcard pattern for input files"
    )
    ccx_solve_parser.add_argument(
        "-p", "--nproc", type=int, default=2, help="Number of processes"
    )
    ccx_solve_parser.add_argument(
        "-c", "--ccxexe", default="ccx", help="Calculix executable"
    )
    ccx_solve_parser.add_argument(
        "-m", "--merged-plies", action="store_true", help="Only process merged plies"
    )

    ccx_post_parser = ccx_subparsers.add_parser("post", help="Postprocess CCX results")
    ccx_post_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )
    ccx_post_parser.add_argument(
        "-w", "--wildcard", default="", help="Wildcard pattern for results"
    )
    ccx_post_parser.add_argument(
        "-n", "--nbins", type=int, default=60, help="Number of bins for tabulation"
    )

    ccx_plot_parser = ccx_subparsers.add_parser("plot", help="Plot CCX results")
    ccx_plot_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )
    # # ccx_plot_parser.add_argument(
    # #     "-w", "--wildcard", default="", help="Wildcard pattern for results"
    # )
    ccx_plot_parser.add_argument(
        "-3",
        "--no-plot3d",
        action="store_false",
        dest="plot3d",
        help="Disable 3D plots",
    )
    ccx_plot_parser.add_argument(
        "-2",
        "--no-plot2d",
        action="store_false",
        dest="plot2d",
        help="Disable 2D plots",
    )

    # 2D subcommands
    twod_parser = subparsers.add_parser("2d", help="2D mesh and ANBA4 operations")
    twod_parser.add_argument("yml", type=Path, help="Path to YAML config file")
    twod_parser.add_argument(
        "-r", "--rotz", type=float, default=0.0, help="Rotation around Z-axis (degrees)"
    )
    twod_parser.add_argument(
        "-P",
        "--no-parallel",
        action="store_false",
        dest="parallel",
        help="Disable parallel processing for mesh2d",
    )
    twod_parser.add_argument(
        "-e", "--anba-env", default="anba4-env", help="Conda environment for ANBA4"
    )
    twod_parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="Output directory for VTU files",
        default=None,
    )
    twod_subparsers = twod_parser.add_subparsers(dest="subcommand", required=False)

    twod_mesh2d_parser = twod_subparsers.add_parser("mesh2d", help="Create 2D meshes")
    twod_mesh2d_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )
    twod_mesh2d_parser.add_argument(
        "-r",
        "--rotz_sub",
        type=float,
        default=0.0,
        dest="rotz",
        help="Rotation around Z-axis (degrees)",
    )
    twod_mesh2d_parser.add_argument(
        "-P",
        "--no-parallel-sub",
        action="store_false",
        dest="parallel",
        help="Disable parallel processing",
    )

    twod_run_anba4_parser = twod_subparsers.add_parser(
        "run-anba4", help="Run ANBA4 on 2D meshes"
    )
    twod_run_anba4_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )
    twod_run_anba4_parser.add_argument(
        "-e",
        "--anba-env-sub",
        default="anba4-env",
        dest="anba_env",
        help="Conda environment for ANBA4",
    )

    twod_clean_parser = twod_subparsers.add_parser("clean", help="Remove msec* files")
    twod_clean_parser.add_argument(
        "yml_sub", type=Path, help="Path to YAML config file", nargs="?", default=None
    )

    # CCBlade subcommand
    ccblade_parser = subparsers.add_parser("ccblade", help="Run CCBlade analysis")
    ccblade_parser.add_argument("yml", type=Path, help="Path to YAML config file")

    # Clean subcommand
    clean_parser = subparsers.add_parser("clean", help="Clean working directory")
    clean_parser.add_argument("yml", type=Path, help="Path to YAML config file")

    args = parser.parse_args()

    state = AppState.get_instance()
    ccx = CcxApp(state, args.yml)
    d2d = TwoDApp(state, args.yml)
    ccb = CCBladeApp(state, args.yml)
    clean = CleanApp(state, args.yml)

    if args.command == "build":
        if not hasattr(args, "subcommand") or args.subcommand is None:
            build = BuildApp(state, args.yml)
            build.build(bondline=args.bondline)
        else:
            build = BuildApp(state, args.yml_sub or args.yml)
            if args.subcommand == "geometry":
                build.geometry()
            elif args.subcommand == "mesh":
                build.mesh()
            elif args.subcommand == "drape":
                build.drape(bondline=args.bondline)
            elif args.subcommand == "mass":
                build.mass()
            elif args.subcommand == "apply-loads":
                build.apply_loads()
    elif args.command == "ccx":
        if not hasattr(args, "subcommand") or args.subcommand is None:
            ccx.ccx(bondline=args.bondline)
        elif args.subcommand == "ccx":
            ccx.ccx(bondline=args.bondline)
        elif args.subcommand == "prep":
            ccx.prep(bondline=args.bondline)
        elif args.subcommand == "solve":
            ccx.solve(
                wildcard=args.wildcard,
                nproc=args.nproc,
                ccxexe=args.ccxexe,
                merged_plies=args.merged_plies,
            )
        elif args.subcommand == "post":
            ccx.post(wildcard=args.wildcard, nbins=args.nbins)
        elif args.subcommand == "plot":
            ccx.plot(
                plot3d=args.plot3d,
                plot2d=args.plot2d,
            )
    elif args.command == "2d":
        if not hasattr(args, "subcommand") or args.subcommand is None:
            d2d.mesh2d(rotz=args.rotz, parallel=args.parallel)
            d2d.run_anba4(anba_env=args.anba_env)
        elif args.subcommand == "mesh2d":
            d2d.mesh2d(args.yml_sub or args.yml, rotz=args.rotz, parallel=args.parallel)
        elif args.subcommand == "run-anba4":
            d2d.run_anba4(args.yml_sub or args.yml, anba_env=args.anba_env)
        elif args.subcommand == "clean":
            d2d.clean()
    elif args.command == "ccblade":
        ccb.ccblade()
    elif args.command == "clean":
        clean.clean()


if __name__ == "__main__":
    main()
