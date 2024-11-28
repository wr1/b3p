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
            self.expand_chamfered_cores()
        return self.dct

    def expand_chamfered_cores(self):
        if self.dct:
            self.dct = build_plybook.expand_chamfered_cores(self.dct)

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
        prefix = wd / dct["general"]["prefix"]

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
        self.geometry(yml)
        self.mesh(yml)
        self.drape(yml, bondline=bondline)
        self.mass(yml)

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

state = AppState.get_instance()

# Register Build Commands
build = BuildApp(state)
build_app.command(build.geometry)
build_app.command(build.mesh)
build_app.command(build.drape)
build_app.command(build.mass)
build_app.command(build.build)
build_app.command(build.mesh2d)

app.command(build_app)

# Register Clean Command
clean = CleanApp(state)
app.command(clean.clean)


def main():
    app()


if __name__ == "__main__":
    main()
