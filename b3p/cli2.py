from pathlib import Path
import os
import glob
import shutil
import pickle
import multiprocessing
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


class AppState:
    """Singleton to manage application state."""

    _state = None

    def __init__(self):
        self.dct = None

    @classmethod
    def get_instance(cls):
        if cls._state is None:
            cls._state = cls()
        return cls._state

    def load_yaml(self, yml: Path):
        if self.dct is None:
            self.dct = yml_portable.yaml_make_portable(yml)
            self.make_workdir(yml)
            self.expand_chamfered_cores()
        return self.dct

    def expand_chamfered_cores(self):
        if self.dct:
            self.dct = build_plybook.expand_chamfered_cores(self.dct)

    def make_workdir(self, yml: Path):
        if self.dct is None:
            self.load_yaml(yml)
        prefix = os.path.join(
            self.dct["general"]["workdir"], self.dct["general"]["prefix"]
        )
        if not os.path.isdir(prefix):
            os.makedirs(prefix)

    def get_prefix(self):
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
        """Build blade geometry."""
        dct = self.state.load_yaml(yml)
        build_blade_geometry.build_blade_geometry(dct)
        prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])
        yml_portable.save_yaml(f"{prefix}_portable.yml", dct)

    def mesh(self, yml: Path):
        """Mesh blade structure."""
        dct = self.state.load_yaml(yml)
        build_blade_structure.build_blade_structure(dct)

    def drape(self, yml: Path, bondline: bool = False):
        """Drape plies onto mesh."""
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
        """Calculate mass of the blade."""
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
        """Build the blade model: geometry, mesh, drape, and mass."""
        dct = self.state.load_yaml(yml)
        self.geometry(yml)
        self.mesh(yml)
        self.drape(yml, bondline=bondline)
        self.mass(yml)
        add_load_to_mesh.add_load_to_mesh(
            dct,
            f"{self.state.get_prefix()}_joined.vtu",
            f"{self.state.get_prefix()}_loads.png",
        )

    def mesh2d(self, yml: Path, rotz=0.0, parallel=True):
        """Create 2D meshes for calculation of 6x6 matrices."""
        dct = self.state.load_yaml(yml)

        if "mesh2d" not in dct:
            print("** No mesh2d section in yml file")
            return

        if "sections" not in dct["mesh2d"]:
            print("** No sections in mesh2d section in yml file")
            return

        sections = dct["mesh2d"]["sections"]
        prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])
        section_meshes = mesh_2d.cut_blade_parallel(
            f"{prefix}_joined.vtu",
            sections,
            if_bondline=False,
            rotz=rotz,
            var=f"{prefix}.var",
            parallel=parallel,
        )
        anba4_prep.anba4_prep(section_meshes)


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
        """Prepare CCX input files."""
        dct = self.state.load_yaml(yml)
        prefix = self.state.get_prefix()
        # os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])

        available_meshes = glob.glob(f"{prefix}_joined.vtu")
        if kwargs.get("bondline"):
            bondline_meshes = glob.glob(f"{prefix}*_bondline.vtu")
            if bondline_meshes:
                available_meshes = bondline_meshes

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
        merged_plies=True,
    ):
        """Run CCX on all input files in the workdir where the loadcase name matches the loadcase wildcard."""
        dct = self.state.load_yaml(yml)

        prefix = self.state.get_prefix()

        if inpfiles is None:
            inpfiles = glob.glob(f"{prefix}*ccx*{wildcard}*.inp")

        inps = [inp for inp in inpfiles if glob.fnmatch.fnmatch(inp, f"*{wildcard}*")]

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
        p.map(os.system, [f"{ccxexe} {inp.replace('.inp','')}" for inp in inps_to_run])
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
        self.state = state

    def mesh2d(self, yml: Path, rotz=0.0, parallel=True):
        """Create 2D meshes for calculation of 6x6 matrices."""
        dct = self.state.load_yaml(yml)

        if "mesh2d" not in dct:
            print("** No mesh2d section in yml file")
            return

        if "sections" not in dct["mesh2d"]:
            print("** No sections in mesh2d section in yml file")
            return

        sections = dct["mesh2d"]["sections"]
        prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])
        section_meshes = mesh_2d.cut_blade_parallel(
            f"{prefix}_joined.vtu",
            sections,
            if_bondline=False,
            rotz=rotz,
            var=f"{prefix}.var",
            parallel=parallel,
        )
        anba4_prep.anba4_prep(section_meshes)


class CleanApp:
    def __init__(self, state):
        self.state = state

    def clean(self, yml: Path):
        """Clean the working directory."""
        dct = self.state.load_yaml(yml)
        prefix = dct["general"]["workdir"]
        if os.path.isdir(prefix):
            shutil.rmtree(prefix)
            print(f"** Removing workdir {prefix}")
        else:
            print(f"** Workdir {prefix} does not exist")


# Initialize Main App
app = App()
build_app = App(name="build")
ccx_app = App(name="ccx")
twod_app = App(name="2d")

state = AppState.get_instance()

# Register Build Commands
build = BuildApp(state)
build_app.command(build.geometry)
build_app.command(build.mesh)
build_app.command(build.drape)
build_app.command(build.mass)
build_app.command(build.mesh2d)
build_app.default(build.build)


# Register CCX Commands
ccx = CcxApp(state)
ccx_app.command(ccx.prep)
ccx_app.command(ccx.solve)
ccx_app.command(ccx.post)
ccx_app.command(ccx.plot)
ccx_app.default(ccx.ccx)

d2d = TwoDApp(state)
twod_app.default(d2d.mesh2d)


# Register Clean Command
clean = CleanApp(state)
app.command(clean.clean)
app.command(build_app)
app.command(ccx_app)
app.command(twod_app)


def main():
    app()


if __name__ == "__main__":
    main()
