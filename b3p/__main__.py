import argparse
from pathlib import Path
from b3p.cli.app_state import AppState
from b3p.cli.build_app import BuildApp
from b3p.cli.ccx_app import CcxApp
from b3p.cli.two_d_app import TwoDApp
from b3p.cli.ccblade_app import CCBladeApp
from b3p.cli.clean_app import CleanApp


def main():
    parser = argparse.ArgumentParser(description="Blade Design CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Build command (top-level)
    build_parser = subparsers.add_parser("build", help="Build the full blade model")
    build_parser.add_argument("yml", type=Path, help="Path to YAML config file")
    build_parser.add_argument(
        "--no-bondline", action="store_false", dest="bondline", help="Exclude bondline"
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
        "--bondline", action="store_true", help="Add bondline to mesh"
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
    ccx_subparsers = ccx_parser.add_subparsers(dest="subcommand", required=True)

    ccx_ccx_parser = ccx_subparsers.add_parser("ccx", help="Run full Calculix process")
    ccx_ccx_parser.add_argument("yml", type=Path, help="Path to YAML config file")
    ccx_ccx_parser.add_argument(
        "--bondline", action="store_true", help="Use bondline meshes"
    )

    ccx_prep_parser = ccx_subparsers.add_parser("prep", help="Prepare CCX input files")
    ccx_prep_parser.add_argument("yml", type=Path, help="Path to YAML config file")
    ccx_prep_parser.add_argument(
        "--bondline", action="store_true", help="Use bondline meshes"
    )

    ccx_solve_parser = ccx_subparsers.add_parser("solve", help="Solve CCX problem")
    ccx_solve_parser.add_argument("yml", type=Path, help="Path to YAML config file")
    ccx_solve_parser.add_argument(
        "--wildcard", default="", help="Wildcard pattern for input files"
    )
    ccx_solve_parser.add_argument(
        "--nproc", type=int, default=2, help="Number of processes"
    )
    ccx_solve_parser.add_argument("--ccxexe", default="ccx", help="Calculix executable")
    ccx_solve_parser.add_argument(
        "--merged-plies", action="store_true", help="Only process merged plies"
    )

    ccx_post_parser = ccx_subparsers.add_parser("post", help="Postprocess CCX results")
    ccx_post_parser.add_argument("yml", type=Path, help="Path to YAML config file")
    ccx_post_parser.add_argument(
        "--wildcard", default="", help="Wildcard pattern for results"
    )
    ccx_post_parser.add_argument(
        "--nbins", type=int, default=60, help="Number of bins for tabulation"
    )

    ccx_plot_parser = ccx_subparsers.add_parser("plot", help="Plot CCX results")
    ccx_plot_parser.add_argument("yml", type=Path, help="Path to YAML config file")
    ccx_plot_parser.add_argument(
        "--wildcard", default="", help="Wildcard pattern for results"
    )
    ccx_plot_parser.add_argument(
        "--no-plot3d", action="store_false", dest="plot3d", help="Disable 3D plots"
    )
    ccx_plot_parser.add_argument(
        "--no-plot2d", action="store_false", dest="plot2d", help="Disable 2D plots"
    )

    # 2D subcommands
    twod_parser = subparsers.add_parser("2d", help="2D mesh and ANBA4 operations")
    twod_subparsers = twod_parser.add_subparsers(dest="subcommand", required=True)

    twod_mesh2d_parser = twod_subparsers.add_parser("mesh2d", help="Create 2D meshes")
    twod_mesh2d_parser.add_argument("yml", type=Path, help="Path to YAML config file")
    twod_mesh2d_parser.add_argument(
        "--rotz", type=float, default=0.0, help="Rotation around Z-axis (degrees)"
    )
    twod_mesh2d_parser.add_argument(
        "--no-parallel",
        action="store_false",
        dest="parallel",
        help="Disable parallel processing",
    )

    twod_run_anba4_parser = twod_subparsers.add_parser(
        "run-anba4", help="Run ANBA4 on 2D meshes"
    )
    twod_run_anba4_parser.add_argument(
        "yml", type=Path, help="Path to YAML config file"
    )
    twod_run_anba4_parser.add_argument(
        "--anba-env", default="anba4-env", help="Conda environment for ANBA4"
    )

    twod_clean_parser = twod_subparsers.add_parser("clean", help="Remove msec* files")
    twod_clean_parser.add_argument("yml", type=Path, help="Path to YAML config file")

    # CCBlade subcommand
    ccblade_parser = subparsers.add_parser("ccblade", help="Run CCBlade analysis")
    ccblade_parser.add_argument("yml", type=Path, help="Path to YAML config file")

    # Clean subcommand
    clean_parser = subparsers.add_parser("clean", help="Clean working directory")
    clean_parser.add_argument("yml", type=Path, help="Path to YAML config file")

    args = parser.parse_args()

    state = AppState.get_instance()
    build = BuildApp(state)
    ccx = CcxApp(state)
    d2d = TwoDApp(state)
    ccb = CCBladeApp(state)
    clean = CleanApp(state)

    if args.command == "build":
        if not hasattr(args, "subcommand") or args.subcommand is None:
            build.build(args.yml, bondline=args.bondline)
        elif args.subcommand == "geometry":
            build.geometry(args.yml_sub or args.yml)
        elif args.subcommand == "mesh":
            build.mesh(args.yml_sub or args.yml)
        elif args.subcommand == "drape":
            build.drape(args.yml_sub or args.yml, bondline=args.bondline)
        elif args.subcommand == "mass":
            build.mass(args.yml_sub or args.yml)
        elif args.subcommand == "apply-loads":
            build.apply_loads(args.yml_sub or args.yml)
    elif args.command == "ccx":
        if args.subcommand == "ccx":
            ccx.ccx(args.yml, bondline=args.bondline)
        elif args.subcommand == "prep":
            ccx.prep(args.yml, bondline=args.bondline)
        elif args.subcommand == "solve":
            ccx.solve(
                args.yml,
                wildcard=args.wildcard,
                nproc=args.nproc,
                ccxexe=args.ccxexe,
                merged_plies=args.merged_plies,
            )
        elif args.subcommand == "post":
            ccx.post(args.yml, wildcard=args.wildcard, nbins=args.nbins)
        elif args.subcommand == "plot":
            ccx.plot(
                args.yml, wildcard=args.wildcard, plot3d=args.plot3d, plot2d=args.plot2d
            )
    elif args.command == "2d":
        if args.subcommand == "mesh2d":
            d2d.mesh2d(args.yml, rotz=args.rotz, parallel=args.parallel)
        elif args.subcommand == "run-anba4":
            d2d.run_anba4(args.yml, anba_env=args.anba_env)
        elif args.subcommand == "clean":
            d2d.clean(args.yml)
    elif args.command == "ccblade":
        ccb.ccblade(args.yml)
    elif args.command == "clean":
        clean.clean(args.yml)


if __name__ == "__main__":
    main()
