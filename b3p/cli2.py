from pathlib import Path

import os
import glob
import shutil
import pickle
import multiprocessing
import subprocess
from cyclopts import App
from b3p import (
    build_blade_geometry,
    build_blade_structure,
    build_plybook,
    drape_mesh,
    combine_meshes,
    mesh_2d,
    add_load_to_mesh,
    mesh2ccx,
    ccx2vtu,
    ccxpost,
    drape_summary,
    anba4_prep,
    add_te_solids,
    yml_portable,
)

try:
    from b3p import ccblade_run

    has_ccblade = True
except ImportError:
    print("** Could not import ccblade_run. Functionality will be disabled.")
    has_ccblade = False


def wslpath_convert(path: str, to_windows: bool = False) -> str:
    """
    Converts paths between WSL and Windows formats using wslpath.

    :param path: The path to convert.
    :param to_windows: Set to True to convert WSL to Windows; False for Windows to WSL.
    :return: The converted path as a string.
    """
    command = ["wslpath", "-w" if to_windows else "-u", path]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    return result.stdout.strip()


class AppState:
    """Singleton to manage application state."""

    _state = None

    def __init__(self):
        self.dct = None

    @classmethod
    def get_instance(cls):
        """
        Returns a singleton instance of the class.

        If the class instance does not already exist, it creates one and stores it
        in the class attribute `_state`. If the instance already exists, it returns
        the existing instance.

        Returns:
            cls: The singleton instance of the class.
        """
        if cls._state is None:
            cls._state = cls()
        return cls._state

    def load_yaml(self, yml: Path):
        """
        Loads and processes a YAML file.

        Args:
            yml (Path): The path to the YAML file to be loaded.

        Returns:
            dict: A dictionary representation of the processed YAML file.
        """
        if self.dct is None:
            self.dct = yml_portable.yaml_make_portable(yml)
            self.make_workdir(yml)
            self.expand_chamfered_cores()
        return self.dct

    def expand_chamfered_cores(self):
        """
        Expands the chamfered cores in the design dictionary.

        This method modifies the `dct` attribute by expanding the chamfered cores
        using the `build_plybook.expand_chamfered_cores` function.

        Returns:
            None
        """
        if self.dct:
            self.dct = build_plybook.expand_chamfered_cores(self.dct)

    def make_workdir(self, yml: Path):
        """
        Creates a working directory based on the provided YAML configuration.

        This method constructs a directory path using the 'workdir' and 'prefix'
        values from the loaded YAML configuration. If the directory does not
        already exist, it will be created.

        Args:
            yml (Path): The path to the YAML file containing the configuration.

        Raises:
            FileNotFoundError: If the YAML file does not exist.
            KeyError: If the required keys ('general', 'workdir', 'prefix') are
                      not present in the YAML configuration.
        """
        if self.dct is None:
            self.load_yaml(yml)
        prefix = os.path.join(
            self.dct["general"]["workdir"], self.dct["general"]["prefix"]
        )
        if not os.path.isdir(prefix):
            os.makedirs(prefix)

    def get_prefix(self):
        """
        Retrieves the prefix path from the configuration dictionary.

        Returns:
            Path or None: The prefix path if the configuration dictionary is not None,
                          otherwise None.
        """
        if self.dct is None:
            return None
        wd = Path(self.dct["general"]["workdir"])
        prefix = wd / self.dct["general"]["prefix"]
        return prefix

    def reset(self):
        self.dct = None


class BuildApp:
    def __init__(self, state):
        self.state = state

    def geometry(self, yml: Path):
        """
        Build blade geometry from a YAML configuration file.

        Args:
            yml (Path): Path to the YAML file containing the blade geometry configuration.

        This method performs the following steps:
        1. Loads the YAML configuration file into a dictionary.
        2. Builds the blade geometry using the configuration dictionary.
        3. Constructs a file path prefix using the 'workdir' and 'prefix' values from the configuration.
        4. Saves the modified configuration dictionary to a new YAML file with a '_portable' suffix.

        Raises:
            FileNotFoundError: If the specified YAML file does not exist.
            KeyError: If required keys are missing in the YAML configuration.
        """
        dct = self.state.load_yaml(yml)
        build_blade_geometry.build_blade_geometry(dct)
        prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])
        yml_portable.save_yaml(f"{prefix}_portable.yml", dct)

    def mesh(self, yml: Path):
        """
        Mesh blade structure.

        Parameters:
        yml (Path): The path to the YAML file containing the blade structure configuration.

        Returns:
        None
        """
        dct = self.state.load_yaml(yml)
        build_blade_structure.build_blade_structure(dct)

    def drape(self, yml: Path, bondline: bool = False):
        """
        Drape plies onto mesh.

        Args:
            yml (Path): Path to the YAML file containing configuration data.
            bondline (bool, optional): If True, add bondline to the mesh. Defaults to False.

        Raises:
            FileNotFoundError: If the plybook file is not found.

        """
        dct = self.state.load_yaml(yml)
        plybookname = "__plybook.pck"
        prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])

        build_plybook.lamplan2plies(dct, plybookname)
        slb = dct["laminates"]["slabs"]
        used_grids = {slb[i]["grid"] for i in slb}

        pbook = os.path.join(dct["general"]["workdir"], plybookname)
        if os.path.exists(pbook):
            plybook = pickle.load(open(pbook, "rb"))
            meshes = []
            for grid in used_grids:
                out = f"{prefix}_{grid}_dr.vtu"
                drape_mesh.drape_mesh(f"{prefix}_{grid}.vtp", plybook, grid, out)
                meshes.append(out)
            combine_meshes.combine_meshes(meshes, f"{prefix}_joined.vtu")
            if bondline:
                add_te_solids.add_bondline(dct)
        else:
            raise FileNotFoundError("Plybook not found")

    def mass(self, yml: Path):
        """
        Calculate the mass of the blade and generate summary tables.

        This method loads a YAML configuration file, checks for the existence of a working directory,
        and if not found, builds a new one. It then calculates the mass of the blade using the drape
        summary and generates CSV and LaTeX files with the mass data.

        Args:
            yml (Path): The path to the YAML configuration file.

        Raises:
            FileNotFoundError: If the working directory specified in the YAML file does not exist.

        Returns:
            None
        """
        dct = self.state.load_yaml(yml)
        wd = Path(dct["general"]["workdir"])
        prefix = self.state.get_prefix()
        if not wd.is_dir():
            print("\n** No workdir found, building new\n")
            self.build(yml)

        mass_table = drape_summary.drape_summary(f"{prefix}_joined.vtu")
        mass_table.to_csv(f"{prefix}_mass.csv")
        mass_table.replace(to_replace="_", value="", regex=True).to_latex(
            f"{prefix}_mass.tex", index=False
        )
        print("Mass table per material")
        print(mass_table)

    def build(self, yml: Path, bondline: bool = True):
        """
        Build the blade model: geometry, mesh, drape, and mass.

        Args:
            yml (Path): The path to the YAML configuration file.
            bondline (bool, optional): A flag indicating whether to include bondline in the drape. Defaults to True.

        Returns:
            None
        """
        dct = self.state.load_yaml(yml)
        self.geometry(yml)
        self.mesh(yml)
        self.drape(yml, bondline=bondline)
        self.mass(yml)
        self.apply_loads(yml)
        # add_load_to_mesh.add_load_to_mesh(
        #     dct,
        #     f"{self.state.get_prefix()}_joined.vtu",
        #     f"{self.state.get_prefix()}_loads.png",
        # )

    def apply_loads(self, yml: Path):
        """
        Reapply loads to the mesh.

        Args:
            yml (Path): Path to the YAML configuration file.

        This method loads the YAML configuration file, retrieves the prefix path, and re-applies the loads
        to the mesh using the 'add_load_to_mesh' function.

        Returns:
            None
        """
        dct = self.state.load_yaml(yml)
        add_load_to_mesh.add_load_to_mesh(
            dct,
            f"{self.state.get_prefix()}_joined.vtu",
            f"{self.state.get_prefix()}_loads.png",
        )


class CcxApp:
    def __init__(self, state):
        self.state = state

    def ccx(self, yml: Path, **kwargs):
        """Run Calculix on this model

        Args:
            yml (Path): Path to the input file.
        """
        self.prep(yml, **kwargs)
        self.solve(yml, **kwargs)
        self.post(yml, **kwargs)
        self.plot(yml, **kwargs)

    def prep(self, yml: Path, **kwargs):
        """
        Prepare CCX input files.

        This method loads a YAML configuration file, determines the appropriate mesh files,
        and generates CCX input files using the mesh2ccx utility.

        Args:
            yml (Path): Path to the YAML configuration file.
            **kwargs: Additional keyword arguments to pass to the mesh2ccx function.

        Keyword Args:
            bondline (bool): If True, use bondline meshes instead of joined meshes.

        Returns:
            None

        Prints:
            A message indicating the output files that were written.
        """
        dct = self.state.load_yaml(yml)
        prefix = self.state.get_prefix()
        # os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])

        available_meshes = glob.glob(f"{prefix}_joined.vtu")
        if kwargs.get("bondline"):
            bondline_meshes = glob.glob(f"{prefix}*_bondline.vtu")
            if bondline_meshes:
                available_meshes = bondline_meshes

        print(f"Available meshes: {available_meshes}")

        output_files = mesh2ccx.mesh2ccx(
            available_meshes[-1],
            matmap=os.path.join(dct["general"]["workdir"], "material_map.json"),
            out=f"{prefix}_ccx.inp",
            **kwargs,
        )
        print(f"Written: {', '.join(output_files)}")

    def solve(
        self,
        yml: Path,
        wildcard="",
        nproc=2,
        ccxexe="ccx",
        inpfiles=None,
        merged_plies=False,
    ):
        """
        Solves the problem defined by the input files using the specified solver.
        Parameters:
        yml (Path): Path to the YAML file containing the state configuration.
        wildcard (str, optional): Wildcard pattern to filter input files. Defaults to "".
        nproc (int, optional): Number of processes to use for parallel execution. Defaults to 2.
        ccxexe (str, optional): Name of the solver executable. Defaults to "ccx".
        inpfiles (list, optional): List of input files to process. If None, files are globbed based on the prefix and wildcard. Defaults to None.
        merged_plies (bool, optional): If True, only process input files with "_mp_" in their name. Defaults to True.
        Returns:
        None
        """

        dct = self.state.load_yaml(yml)

        prefix = self.state.get_prefix()

        if inpfiles is None:
            inpfiles = glob.glob(f"{prefix}*ccx*{wildcard}*.inp")

        inps = [inp for inp in inpfiles]
        # if glob.fnmatch.fnmatch(inp, f"*{wildcard}*")]

        if merged_plies:
            inps = [inp for inp in inps if "_mp_" in inp]

        if inps == []:
            print(f"** No inps found matching {prefix}*{wildcard}*inp")
            return

        inps_to_run = []
        for inp in inps:
            frd_file = inp.replace(".inp", ".frd")
            # check if frd exists
            if not os.path.exists(frd_file):
                inps_to_run.append(inp)
            else:
                # check if the frd is a complete file
                with open(frd_file, "rb") as f:
                    f.seek(-5, 2)
                    y = f.read()
                    if y != b"9999\n":
                        inps_to_run.append(inp)
                    else:
                        print(f"** Skipping {inp} because {frd_file} exists")

        print(f"running ccx on: {inps_to_run} using {wildcard}")
        p = multiprocessing.Pool(nproc)
        p.map(os.system, [f"{ccxexe} {inp.replace('.inp', '')}" for inp in inps_to_run])
        p.close()

    def post(self, yml: Path, wildcard="", nbins=60):
        """Postprocess CCX results."""
        dct = self.state.load_yaml(yml)
        ccxpost = ccx2vtu.ccx2vtu(dct["general"]["workdir"], wildcard=wildcard)
        ccxpost.load_grids()
        ccxpost.tabulate(nbins)

    def plot(self, yml: Path, wildcard="", plot3d=True, plot2d=True):
        """Plot CCX results."""
        dct = self.state.load_yaml(yml)
        plotter = ccxpost.plot_ccx(dct["general"]["workdir"], wildcard=wildcard)
        if plot3d:
            plotter.plot3d()
        if plot2d:
            plotter.plot2d(wildcard=wildcard)


class TwoDApp:
    def __init__(self, state):
        """
        Initialize the instance with the given state.
        Args:
            state: The initial state to be assigned to the instance.
        """
        self.state = state

    def mesh2d(self, yml: Path, rotz=0.0, parallel=True):
        """
        Create 2D meshes for calculation of 6x6 matrices.

        Args:
            yml (Path): Path to the YAML file containing mesh configuration.
            rotz (float, optional): Rotation around the Z-axis in degrees. Defaults to 0.0.
            parallel (bool, optional): Flag to enable parallel processing. Defaults to True.

        Returns:
            list: A list of section meshes prepared for anba4 calculations.

        Raises:
            FileNotFoundError: If the YAML file does not exist.
            KeyError: If required sections are missing in the YAML file.
        """
        dct = self.state.load_yaml(yml)

        if "mesh2d" not in dct:
            print("** No mesh2d section in yml file")
            return

        if "sections" not in dct["mesh2d"]:
            print("** No sections in mesh2d section in yml file")
            return

        sections = dct["mesh2d"]["sections"]
        yml_dir = yml.parent
        prefix = os.path.join(
            yml_dir, dct["general"]["workdir"], dct["general"]["prefix"]
        )
        # prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])
        section_meshes = mesh_2d.cut_blade_parallel(
            f"{prefix}_joined.vtu",
            sections,
            if_bondline=False,
            rotz=rotz,
            var=f"{prefix}.var",
            parallel=parallel,
        )
        return anba4_prep.anba4_prep(section_meshes)

    def run_anba4(self, yml: Path, meshes: list = None, anba_env="anba4-env"):
        """
        Run ANBA4 on the meshes (2D or single).

        Parameters:
        yml (Path): Path to the YAML configuration file.
        meshes (list, optional): List of mesh file paths. If None, meshes are extracted from the YAML.
        anba_env (str): Name of the conda environment to use for running ANBA4. Default is "anba4-env".

        Returns:
        None

        This method performs the following steps:
        1. Checks if Conda is installed on the system.
        2. Verifies if the specified Conda environment exists.
        3. Loads the YAML configuration file.
        4. Extracts the meshes from the YAML configuration file if not provided.
        5. Constructs the path to the material map JSON file.
        6. Runs the ANBA4 solver using the specified Conda environment and the extracted meshes and material map.
        """
        # check if system has conda
        conda_path = os.environ.get("CONDA_EXE") or shutil.which("conda")
        if conda_path is None:
            print(
                "** Conda not found - please install conda.  Ensure CONDA_EXE is set if Conda is not in your PATH."
            )
            return

        # check if the conda environment exists
        result = subprocess.run(
            [conda_path, "env", "list"], capture_output=True, text=True
        )
        if result.returncode != 0 or anba_env not in result.stdout:
            print(f"** Conda environment {anba_env} not found - please create it")
            return
        else:
            print(f"** Using Conda environment for running anba4 {anba_env}")

        dct = self.state.load_yaml(yml)

        if meshes is None:
            meshes = self.mesh2d(yml)

        yml_dir = yml.parent
        material_map = os.path.join(
            yml_dir, dct["general"]["workdir"], "material_map.json"
        )

        # Call the new CLI script using the conda environment
        script_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "anba4_solve.py")
        )
        print(f"** conda path: {conda_path}")
        print(f"** Running ANBA4 using {script_path}")

        return subprocess.run(
            [
                conda_path,
                "run",
                "-n",
                anba_env,
                "python",
                script_path,
                *meshes,
                material_map,
            ],
            env={
                **os.environ.copy(),
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "OMP_NUM_THREADS": "1",
                "CUDA_VISIBLE_DEVICES": "-1",
            },
        ).returncode

    def clean(self, yml: Path):
        """
        Remove all msec* files in the working directory.

        Args:
            yml (Path): Path to the YAML configuration file.

        This method loads the YAML file to retrieve the working directory path.
        It then removes all files in the directory that start with 'msec'.
        """
        dct = self.state.load_yaml(yml)
        workdir = Path(dct["general"]["workdir"])
        for msec_file in workdir.glob("msec*"):
            try:
                msec_file.unlink()
                print(f"Removed {msec_file}")
            except Exception as e:
                print(f"Failed to remove {msec_file}: {e}")


class CCBladeApp:
    def __init__(self, state):
        self.state = state

    def ccblade(self, yml: Path):
        """
        Run CCBlade on this model.

        Args:
            yml (Path): Path to the input file.
        """
        if has_ccblade:
            # Initialize ccblade_run with the blade file
            ccblade = ccblade_run.ccblade_run(yml)

            # Run the CCBlade analysis
            ccblade.run()
        else:
            print("** ccblade_run is not available.")


class CleanApp:
    def __init__(self, state):
        """
        Initializes a new instance of the class with the given state.

        Args:
            state: The initial state to be assigned to the instance.
        """
        self.state = state

    def clean(self, yml: Path):
        """
        Clean the working directory specified in the YAML configuration file.

        Args:
            yml (Path): Path to the YAML configuration file.

        The method loads the YAML file to retrieve the working directory path.
        If the directory exists, it is removed along with all its contents.
        If the directory does not exist, a message is printed indicating so.
        """
        dct = self.state.load_yaml(yml)
        prefix = dct["general"]["workdir"]
        if os.path.isdir(prefix):
            shutil.rmtree(prefix)
            print(f"** Removing workdir {prefix}")
        else:
            print(f"** Workdir {prefix} does not exist")


# Initialize Main App
app = App(help="""Blade Design CLI""")
build_app = App(name="build")
ccx_app = App(name="ccx")
twod_app = App(name="2d")
ccblade_app = App(name="ccblade")

state = AppState.get_instance()

# Register Build Commands
build = BuildApp(state)
build_app.command(build.geometry)
build_app.command(build.mesh)
build_app.command(build.drape)
build_app.command(build.mass)
build_app.command(build.apply_loads)
build_app.default(build.build)


# Register CCX Commands
ccx = CcxApp(state)
ccx_app.command(ccx.prep)
ccx_app.command(ccx.solve)
ccx_app.command(ccx.post)
ccx_app.command(ccx.plot)
ccx_app.default(ccx.ccx)

d2d = TwoDApp(state)
twod_app.command(d2d.mesh2d)
twod_app.default(d2d.run_anba4)
twod_app.command(d2d.clean)

ccb = CCBladeApp(state)
ccblade_app.default(ccb.ccblade)


# Register Clean Command
clean = CleanApp(state)
app.command(clean.clean)
app.command(build_app)
app.command(ccx_app)
app.command(twod_app)
app.command(ccblade_app)


def main():
    app()


if __name__ == "__main__":
    main()
