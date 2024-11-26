#! /usr/bin/env python

import cyclopts
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
import os
import pickle
import shutil
import glob
import multiprocessing
from pathlib import Path

app = cyclopts.App()

buildapp = cyclopts.App(name="build")
m2dapp = cyclopts.App(name="2d")
ccxapp = cyclopts.App(name="ccx")


app.command(buildapp)
app.command(m2dapp)
app.command(ccxapp)


dct = None


def get_prefix():
    global dct
    return os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])


def has_mesh(dct):
    prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])
    return glob.glob(f"{prefix}*.vtu") != []


@app.command
def clean(yml: Path):
    """Clean the working directory.

    Args:
        yml (Path): b3p input file
    """
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    prefix = dct["general"]["workdir"]  # , dct["general"]["prefix"])
    if os.path.isdir(prefix):
        shutil.rmtree(prefix)
        print(f"** Removing workdir {prefix}")
    else:
        print(f"** Workdir {prefix} does not exist")


@buildapp.command
def geometry(yml: Path):
    """Build blade geometry

    Args:
        yml (Path): b3p input file
    """

    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])
    build_blade_geometry.build_blade_geometry(dct)
    yml_portable.save_yaml(f"{prefix}_portable.yml", dct)


@buildapp.command
def mesh(yml: Path):
    """Mesh blade

    Args:
        yml (Path): b3p input file
    """
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    build_blade_structure.build_blade_structure(dct)


@buildapp.command
def drape(yml: Path, bondline: bool = False):
    """Drape plies onto mesh

    Args:
        yml (Path): b3p input file
    """
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

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


@buildapp.command
def mass(yml: Path):
    """Calculate mass of the blade, requires drape."""
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    wd = Path(dct["general"]["workdir"])
    prefix = wd / dct["general"]["prefix"]

    if not wd.is_dir():
        # exit("** Workdir does not exist, cannot compute mass")
        print("\n" * 10 + "** no workdir found, building new" + "\n" * 10)
        build(yml)

    mass_table = drape_summary.drape_summary(f"{prefix}_joined.vtu")
    mass_table.to_csv(f"{prefix}_mass.csv")
    mass_table.replace(to_replace="_", value="", regex=True).to_latex(
        f"{prefix}_mass.tex", index=False
    )
    print("Mass table per material")
    print(mass_table)


def getdct(yml: Path):
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)
    return dct
    # return has_mesh(dct)


@buildapp.default
def build(yml: Path, bondline: bool = True):
    """Build the blade model, geometry, mesh, drape

    Args:
        yml (Path): b3p input file
    """
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    geometry(yml)
    mesh(yml)
    drape(yml, bondline=bondline)
    mass(yml)
    prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])
    add_load_to_mesh.add_load_to_mesh(
        dct, f"{prefix}_joined.vtu", f"{prefix}_loads.png"
    )


@m2dapp.command
def mesh2d(yml: Path, rotz=0.0, parallel=True):
    """Create 2d meshes for calculation of 6x6 matrices

    Args:
        yml (Path): b3p input file
        rotz (float, optional): Rotation angle. Defaults to 0.0.
        parallel (bool, optional): Use parallel processing. Defaults to True.
    """
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    if not has_mesh(dct):
        build(yml)

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


@m2dapp.default
def m2d(yml: Path):
    """Create 2d meshes for calculation of 6x6 matrices

    Args:
        yml (Path): b3p input file
    """
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    mesh2d(yml)


@ccxapp.command
def ccxprep(
    yml: Path,
    merge_adjacent_layers=True,
    single_step=False,
    quadratic=True,
    force_isotropic=False,
    zeroangle=False,
    add_centers=False,
    export_plygroups=False,
    export_hyperworks=False,
    meshonly=False,
    buckling=False,
    bondline=False,
):
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])

    if not has_mesh(dct):
        build(yml)

    if bondline:
        available_meshes = glob.glob(f"{prefix}*_bondline.vtu")

    if not bondline or available_meshes == []:
        available_meshes = glob.glob(f"{prefix}*_joined.vtu")

    output_files = mesh2ccx.mesh2ccx(
        available_meshes[-1],
        matmap=os.path.join(dct["general"]["workdir"], "material_map.json"),
        out=f"{prefix}_ccx.inp",
        merge_adjacent_layers=merge_adjacent_layers,
        single_step=single_step,
        quadratic=quadratic,
        force_isotropic=force_isotropic,
        zeroangle=zeroangle,
        add_centers=add_centers,
        export_plygroups=export_plygroups,
        buckling=buckling,
        meshonly=meshonly,
        export_hyperworks=export_hyperworks,
    )
    print(f"written " + "\n".join(output_files))
    return output_files


@ccxapp.command
def ccxsolve(yml: Path, wildcard="", nproc=2, ccxexe="ccx", inpfiles=[]):
    """Run ccx on all inp files in workdir that match wildcard

    Args:
        yml (Path): b3p input file
        wildcard (str, optional): wildcard for inp files. Defaults to "".
        nproc (int, optional): number of processors. Defaults to 2.
        ccxexe (str, optional): path to ccx executable. Defaults to "ccx".
        inpfiles (list, optional): list of inp files to run. Defaults to [].
    """
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    prefix = os.path.join(dct["general"]["workdir"], dct["general"]["prefix"])
    if inpfiles == []:
        inps = glob.glob(f"{prefix}*{wildcard}*inp")
    else:
        inps = inpfiles

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


@ccxapp.command
def ccx_post(yml: Path, wildcard="", nbins=60):
    """Postprocess ccx results

    Args:
        yml (Path): b3p input file
        wildcard (str, optional): wildcard for inp files. Defaults to "".
        nbins (int, optional): number of bins. Defaults to 60.
    """
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    ccxpost = ccx2vtu.ccx2vtu(dct["general"]["workdir"], wildcard=wildcard)
    ccxpost.load_grids()
    ccxpost.tabulate(nbins)


@ccxapp.command
def ccxplot(yml: Path, wildcard="", plot3d=True, plot2d=True):
    """Plot ccx results

    Args:
        yml (Path): b3p input file
        wildcard (str, optional): wildcard for inp files. Defaults to "".
        plot3d (bool, optional): plot 3d. Defaults to True.
        plot2d (bool, optional): plot 2d. Defaults to True.
    """
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    plotter = ccxpost.plot_ccx(dct["general"]["workdir"], wildcard=wildcard)
    if plot3d:
        plotter.plot3d()
    if plot2d:
        plotter.plot2d(wildcard=wildcard)


@ccxapp.default
def ccx(yml: Path, wildcard=""):
    """Run the ccx analysis

    Args:
        yml (Path): b3p input file
    """
    global dct
    if dct is None:
        dct = yml_portable.yaml_make_portable(yml)

    ccxprep(yml)
    ccxsolve(yml, wildcard=wildcard)
    ccx_post(yml)
    ccxplot(yml)


@app.command
def shell():
    """
    Start an interactive shell.
    """
    app.interactive_shell(quit=["q", "quit", "exit"])


# class cli:
#     """Fire command line interface for b3p"""

#     def __init__(self, yml: Path ):
#         dct = yml_portable.yaml_make_portable(yml)
#         plybookname = "__plybook.pck"
#         prefix = os.path.join(
#             dct["general"]["workdir"],
#             dct["general"]["prefix"],
#         )

#     def clean(self):
#         """Clean up workdir"""
#         wd = self.dct["general"]["workdir"]
#         if os.path.isdir(wd):
#             shutil.rmtree(wd)
#             print(f"** Removing workdir {wd}")
#         else:
#             print(f"** Workdir {wd} does not exist")
#         return self

#     def __str__(self):
#         return ""

#     def geometry(self):
#         """Build blade geometry based on yaml input file"""
#         build_blade_geometry.build_blade_geometry(self.dct)
#         yml_portable.save_yaml(f"{self.prefix}_portable.yml", self.dct)
#         # yml_portable.save_yaml_portable(f"{self.prefix}_portable.yml")
#         return self

#     def mesh(self):
#         """Build blade mesh based on yaml input file, requires blade geometry"""
#         build_blade_structure.build_blade_structure(self.dct)
#         return self

#     def __plybook(self):
#         """make plybook"""
#         build_plybook.lamplan2plies(self.dct, self.plybookname)
#         return self

#     def drape(self):
#         """Assign plies in plybook to elements in the structural mesh"""
#         self.__plybook()
#         slb = self.dct["laminates"]["slabs"]
#         used_grids = {slb[i]["grid"] for i in slb}

#         pbook = os.path.join(self.dct["general"]["workdir"], self.plybookname)

#         if os.path.exists(pbook):
#             plybook = pickle.load(open(pbook, "rb"))
#             meshes = []
#             for grid in used_grids:
#                 out = f"{self.prefix}_{grid}_dr.vtu"
#                 drape_mesh.drape_mesh(f"{self.prefix}_{grid}.vtp", plybook, grid, out)
#                 meshes.append(out)
#             combine_meshes.combine_meshes(meshes, f"{self.prefix}_joined.vtu")
#         else:
#             raise FileNotFoundError("Plybook not found")
#         return self

#     def mesh2d(self, rotz=0.0, parallel=True):
#         """Create 2d meshes for calculation of 6x6 stiffness and matrices"""
#         if "mesh2d" not in self.dct:
#             print("** No mesh2d section in yml file")
#             return self
#         if "sections" not in self.dct["mesh2d"]:
#             print("** No sections in mesh2d section in yml file")
#             return self

#         sections = self.dct["mesh2d"]["sections"]
#         section_meshes = mesh_2d.cut_blade_parallel(
#             f"{self.prefix}_joined.vtu",
#             sections,
#             if_bondline=False,
#             rotz=rotz,
#             var=f"{self.prefix}.var",
#             parallel=parallel,
#         )
#         anba4_prep.anba4_prep(section_meshes)
#         return self

#     def show(self):
#         return self

#     def bondline(self):
#         add_te_solids.add_bondline(self.dct)

#     def ccxprep(
#         self,
#         merge_adjacent_layers=True,
#         single_step=False,
#         quadratic=True,
#         force_isotropic=False,
#         zeroangle=False,
#         add_centers=False,
#         export_plygroups=False,
#         export_hyperworks=False,
#         meshonly=False,
#         buckling=False,
#         bondline=False,
#     ):
#         print("** create ccx input file")

#         if bondline:
#             available_meshes = glob.glob(f"{self.prefix}*_bondline.vtu")

#         if not bondline or available_meshes == []:
#             available_meshes = glob.glob(f"{self.prefix}*_joined.vtu")

#         # print("\tavailable_meshes", available_meshes)

#         output_files = mesh2ccx.mesh2ccx(
#             available_meshes[-1],
#             matmap=os.path.join(self.dct["general"]["workdir"], "material_map.json"),
#             out=f"{self.prefix}_ccx.inp",
#             merge_adjacent_layers=merge_adjacent_layers,
#             single_step=single_step,
#             quadratic=quadratic,
#             force_isotropic=force_isotropic,
#             zeroangle=zeroangle,
#             add_centers=add_centers,
#             export_plygroups=export_plygroups,
#             buckling=buckling,
#             meshonly=meshonly,
#             export_hyperworks=export_hyperworks,
#         )
#         return output_files

#     def ccxsolve(self, wildcard="", nproc=2, ccxexe="ccx", inpfiles=[]):
#         """Run ccx on all inp files in workdir that match wildcard"""

#         if inpfiles == []:
#             inps = glob.glob(f"{self.prefix}*{wildcard}*inp")
#         else:
#             inps = inpfiles

#         if inps == []:
#             print(f"** No inps found matching {self.prefix}*{wildcard}*inp")
#             return self

#         inps_to_run = []
#         for inp in inps:
#             frd_file = inp.replace(".inp", ".frd")
#             # check if frd exists
#             if not os.path.exists(frd_file):
#                 inps_to_run.append(inp)
#             else:
#                 # check if the frd is a complete file
#                 with open(frd_file, "rb") as f:
#                     f.seek(-5, 2)
#                     y = f.read()
#                     if y != b"9999\n":
#                         inps_to_run.append(inp)
#                     else:
#                         print(f"** Skipping {inp} because {frd_file} exists")

#         print(f"running ccx on: {inps_to_run} using {wildcard}")
#         p = multiprocessing.Pool(nproc)
#         p.map(os.system, [f"{ccxexe} {inp.replace('.inp','')}" for inp in inps_to_run])
#         p.close()
#         return self

#     def ccxpost(self, wildcard="", nbins=60):
#         """Postprocess ccx results"""
#         post = ccx2vtu.ccx2vtu(self.dct["general"]["workdir"], wildcard=wildcard)
#         post.load_grids()
#         post.tabulate(nbins)
#         return self

#     def ccxplot(self, wildcard="", plot3d=True, plot2d=True):
#         plotter = ccxpost.plot_ccx(self.dct["general"]["workdir"], wildcard=wildcard)
#         if plot3d:
#             plotter.plot3d()
#         if plot2d:
#             plotter.plot2d(wildcard=wildcard)
#         return self

#     def build(self):
#         """Build the whole model, geometry, mesh, drape"""
#         self.dct = build_plybook.expand_chamfered_cores(self.dct)
#         self.geometry()
#         self.mesh()
#         self.drape()
#         add_load_to_mesh.add_load_to_mesh(
#             self.dct, f"{self.prefix}_joined.vtu", f"{self.prefix}_loads.png"
#         )
#         return self

#     def mass(self):
#         """Calculate mass of the blade, requires drape."""
#         mass_table = drape_summary.drape_summary(f"{self.prefix}_joined.vtu")
#         mass_table.to_csv(f"{self.prefix}_mass.csv")
#         mass_table.replace(to_replace="_", value="", regex=True).to_latex(
#             f"{self.prefix}_mass.tex", index=False
#         )
#         return self


# # def main():
# #     fire.core.Display = lambda lines, out: print(*lines, file=out)
# #     fire.Fire(cli)


def main():
    app()


if __name__ == "__main__":
    main()
