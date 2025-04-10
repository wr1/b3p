from pathlib import Path
import os
import glob
import multiprocessing
from b3p.ccx import mesh2ccx, ccx2vtu, ccxpost


class CcxApp:
    def __init__(self, state):
        self.state = state

    def ccx(self, yml: Path, **kwargs):
        self.prep(yml, **kwargs)
        self.solve(yml, **kwargs)
        self.post(yml, **kwargs)
        self.plot(yml, **kwargs)

    def prep(self, yml: Path, bondline=False, **kwargs):
        dct = self.state.load_yaml(yml)
        prefix = self.state.get_prefix()
        available_meshes = glob.glob(f"{prefix}_joined.vtu")
        print(f"Available meshes: {prefix}")
        if bondline:
            bondline_meshes = glob.glob(f"{prefix}*_bondline.vtu")
            if bondline_meshes:
                available_meshes = bondline_meshes

        print(f"Available meshes: {available_meshes}")
        if available_meshes == []:
            print("** No meshes found, did you build the blade geometry?")
            return

        output_files = mesh2ccx.mesh2ccx(
            available_meshes[-1],
            matmap=os.path.join(dct["general"]["workdir"], "material_map.json"),
            out=f"{prefix}_ccx.inp",
            bondline=bondline,
            **{k: v for k, v in kwargs.items() if k != "bondline"},  # Pass other kwargs
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
        bondline=False,  # Added to accept bondline, even if unused here
        **kwargs,  # Accept additional kwargs to avoid errors
    ):
        dct = self.state.load_yaml(yml)
        prefix = self.state.get_prefix()
        if inpfiles is None:
            inpfiles = glob.glob(f"{prefix}*ccx*{wildcard}*.inp")
        inps = [inp for inp in inpfiles]
        if merged_plies:
            inps = [inp for inp in inps if "_mp_" in inp]
        if not inps:
            print(f"** No inps found matching {prefix}*{wildcard}*inp")
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
                        print(f"** Skipping {inp} because {frd_file} exists")
        print(f"running ccx on: {inps_to_run} using {wildcard}")
        p = multiprocessing.Pool(nproc)
        p.map(os.system, [f"{ccxexe} {inp.replace('.inp', '')}" for inp in inps_to_run])
        p.close()

    def post(self, yml: Path, wildcard="", nbins=60, bondline=False, **kwargs):
        dct = self.state.load_yaml(yml)
        ccxpost = ccx2vtu.ccx2vtu(dct["general"]["workdir"], wildcard=wildcard)
        ccxpost.load_grids()
        ccxpost.tabulate(nbins)

    def plot(
        self, yml: Path, wildcard="", plot3d=True, plot2d=True, bondline=False, **kwargs
    ):
        dct = self.state.load_yaml(yml)
        plotter = ccxpost.plot_ccx(dct["general"]["workdir"], wildcard=wildcard)
        if plot3d:
            plotter.plot3d()
        if plot2d:
            plotter.plot2d(wildcard=wildcard)
