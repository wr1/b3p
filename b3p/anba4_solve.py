#! /usr/bin/env python3
import pyvista as pv
from dolfin import *
from anba4 import *
from anba4.voight_notation import *
import argparse
import multiprocessing
import json
import numpy as np
import yaml
import os
from functools import partial


def get_material_db(material_map):
    assert os.path.isfile(material_map)
    mm = json.load(open(material_map, "r"))
    gdir = os.path.dirname(material_map)

    mat_db = None
    if "matdb" in mm:  # check if the material map file points to a material db
        mat_db = yaml.safe_load(open(os.path.join(gdir, mm["matdb"])))

        # check if there is a -1 material in the matdb, and assign it
        # this material ID is used by the section mesher for the bondlines
        # that connect the webs to the shell
        if "-1" in mat_db:
            mm["-1"] = -1
    else:
        exit(
            "material map available, but no link to material db, need matdb definition to do FEA"
        )

    mm_inv = {v: k for k, v in mm.items()}

    materials = {}
    for i in mm_inv:
        if i != mm["matdb"]:
            matdb_entry = mat_db[mm_inv[i]]
            if "tEx" in matdb_entry:  # ortho material
                matMechanicProp = np.zeros((3, 3))
                matMechanicProp[0, 0] = matdb_entry["tEz"] * 1e6  # e_xx
                matMechanicProp[0, 1] = matdb_entry["tEy"] * 1e6  # e_yy
                matMechanicProp[0, 2] = matdb_entry["tEx"] * 1e6  # e_zz

                matMechanicProp[1, 0] = matdb_entry["tGxz"] * 1e6  # g_yz
                matMechanicProp[1, 1] = matdb_entry["tGxy"] * 1e6  # g_xz
                matMechanicProp[1, 2] = matdb_entry["tGyz"] * 1e6  # g_xy

                matMechanicProp[2, 0] = matdb_entry["tnuxz"]  # nu_zy
                matMechanicProp[2, 1] = matdb_entry["tnuxy"]  # nu_zx
                matMechanicProp[2, 2] = matdb_entry["tnuyz"]  # nu_xy

                materials[i] = material.OrthotropicMaterial(
                    matMechanicProp, matdb_entry["rho"]
                )
            else:
                materials[i] = material.IsotropicMaterial(
                    [
                        (
                            matdb_entry["E"] * 1e6
                            if "E" in matdb_entry
                            else matdb_entry["Ex"] * 1e6
                        ),
                        matdb_entry["nu"],
                    ],
                    matdb_entry["rho"],
                )

    return materials


def export_unit_strains(anba, result_file_name):
    result_file = XDMFFile(result_file_name)
    result_file.parameters["functions_share_mesh"] = True
    result_file.parameters["rewrite_function_mesh"] = False
    result_file.parameters["flush_output"] = True

    # anba.stress_field([1., 0., 0.,], [0., 0., 0.], "local", "paraview")
    anba.strain_field([0, 0, 0], [1, 0, 0], "local", "paraview")

    result_file.write(anba.STRAIN, t=0.0)

    anba.strain_field([0, 0, 0], [0, 1, 0], "local", "paraview")

    result_file.write(anba.STRAIN, t=1.0)

    anba.strain_field([0, 0, 0], [0, 0, 1], "local", "paraview")

    result_file.write(anba.STRAIN, t=2.0)

    print(f"writing results to {result_file_name}")


def run_mesh(meshname, matdb_name):
    print(f"run {meshname}")
    matdb = get_material_db(matdb_name)

    infile = XDMFFile(meshname)

    mesh = Mesh()
    infile.read(mesh)

    pvmesh = pv.read(meshname)

    # Basic material parameters. 9 is needed for orthotropic materials.
    # TODO materials and orientations
    materials = MeshFunction("size_t", mesh, mesh.topology().dim())
    fiber_orientations = MeshFunction("double", mesh, mesh.topology().dim())
    plane_orientations = MeshFunction("double", mesh, mesh.topology().dim())

    # material indices (from material_map)
    matuniq = np.unique(pvmesh.cell_data["mat"])

    # map for materials in anba (index in matuniq)
    mat_map_0 = dict(zip(matuniq, range(len(matuniq))))

    matids = [mat_map_0[i] for i in pvmesh.cell_data["mat"].tolist()]

    plane_angles = list(pvmesh.cell_data["angle2"])

    materials.set_values(matids)

    # TODO, doesn't work for off axis laminates for now
    fiber_orientations.set_all(0.0)

    # transverse orientations
    plane_orientations.set_values(plane_angles)

    # Build material property library.

    matLibrary = [matdb[i] for i in matuniq]

    anba = anbax(mesh, 2, matLibrary, materials, plane_orientations, fiber_orientations)
    stiff = anba.compute()

    resfilename = meshname.replace(".xdmf", "_results.xdmf")

    export_unit_strains(anba, resfilename)

    stiffness_matrix = stiff.getDenseArray()

    mass = anba.inertia()

    mass_matrix = mass.getDenseArray()

    decoupled_stiff = DecoupleStiffness(stiff)

    angle = PrincipalAxesRotationAngle(decoupled_stiff)

    mass_center = ComputeMassCenter(mass)
    tension_center = ComputeTensionCenter(stiffness_matrix)
    shear_center = ComputeShearCenter(stiffness_matrix)
    output = {
        "name": meshname,
        "stiffness": stiffness_matrix.tolist(),
        "mass_matrix": mass_matrix.tolist(),
        "decoupled_stiffness": decoupled_stiff.tolist(),
        "principal_axes_rotation": angle,
        "mass_center": mass_center,
        "tension_center": tension_center,
        "shear_center": shear_center,
    }

    with open(f"{meshname}.json", "w") as write_file:
        json.dump(output, write_file, indent=4)


def main():
    p = argparse.ArgumentParser(description="run a series of sections through anba4")
    p.add_argument("meshes", nargs="*")
    p.add_argument("matdb", help="material map json file")
    p.add_argument("--debug", action="store_true")
    args = p.parse_args()

    print(args.meshes)

    part = partial(run_mesh, args.matdb)
    if args.debug:
        for i in args.meshes:
            part(i)
    else:
        processes = []
        for mesh in args.meshes:
            p = multiprocessing.Process(target=part, args=(mesh,))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()


if __name__ == "__main__":
    main()
