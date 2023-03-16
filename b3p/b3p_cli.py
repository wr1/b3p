#! /usr/bin/env python

import fire
from b3p import (
    build_blade_geometry,
    build_blade_structure,
    build_plybook,
    drape_mesh,
    combine_meshes,
    mesh_2d,
    add_load_to_mesh,
    mesh2ccx,
)
from utils import yml_portable
import os
import pickle
import numpy as np
import shutil
import scipy as sp
from copy import deepcopy


class cli:
    """Fire command line interface for b3p"""

    def __init__(self, yml):
        self.dct = yml_portable.yaml_make_portable(yml, True)
        self.plybookname = "__plybook.pck"
        self.prefix = os.path.join(
            self.dct["general"]["workdir"], self.dct["general"]["prefix"]
        )

    def clean(self):
        """Clean up workdir"""
        if os.path.isdir(self.dct["general"]["workdir"]):
            shutil.rmtree(self.dct["general"]["workdir"])
            print(f"** Removing workdir {self.dct['general']['workdir']}")
        else:
            print(f"** Workdir {self.dct['general']['workdir']} does not exist")
        return self

    def __str__(self):
        return ""

    def geometry(self):
        """Build blade geometry based on yaml input file"""
        build_blade_geometry.build_blade_geometry(self.dct)
        return self

    def mesh(self):
        """Build blade mesh based on yaml input file, requires blade geometry"""
        build_blade_structure.build_blade_structure(self.dct)
        return self

    def __plybook(self):
        """make plybook"""
        build_plybook.lamplan2plies(self.dct, self.plybookname)
        return self

    def drape(self):
        """Assign plies in plybook to elements in the structural mesh"""
        self.__plybook()
        slb = self.dct["laminates"]["slabs"]
        used_grids = {slb[i]["grid"] for i in slb}

        pbook = os.path.join(self.dct["general"]["workdir"], self.plybookname)

        if os.path.exists(pbook):
            plybook = pickle.load(open(pbook, "rb"))
            meshes = []
            for grid in used_grids:
                out = f"{self.prefix}_{grid}_dr.vtu"
                drape_mesh.drape_mesh(f"{self.prefix}_{grid}.vtp", plybook, grid, out)
                meshes.append(out)
            combine_meshes.combine_meshes(meshes, f"{self.prefix}_joined.vtu")
        else:
            raise FileNotFoundError("Plybook not found")
        return self

    def mesh2d(self, z_start=0, z_end=100, nsec=100):
        """Create 2d meshes for calculation of 6x6 stiffness and matrices"""
        mesh_2d.cut_blade_parallel(
            f"{self.prefix}_joined.vtu",
            np.linspace(z_start, z_end, nsec).tolist(),
            if_bondline=False,
            rotz=0.0,
            var=f"{self.prefix}.var",
        )
        return self

    def show(self):
        print(self)
        return self

    def ccx(self):
        grid = add_load_to_mesh.add_load_to_mesh(self.dct, f"{self.prefix}_joined.vtu")
        mesh2ccx.mesh2ccx(
            f"{self.prefix}_joined.vtu",
            matmap=os.path.join(self.dct["general"]["workdir"], "material_map.json"),
            out=f"{self.prefix}_ccx.inp",
            merge_adjacent_layers=True,
            single_step=True,
            quadratic=True,
        )
        return self

    def all(self):
        """Run all steps, geometry, mesh, drape, mesh2d"""
        self.geometry()
        self.mesh()
        self.drape()
        self.mesh2d(z_start=0.0001, z_end=100, nsec=50)


def main():
    fire.Fire(cli)


if __name__ == "__main__":
    main()
