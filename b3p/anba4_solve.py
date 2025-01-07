#! /usr/bin/env python3
import pyvista as pv
from dolfin import Mesh, XDMFFile, MeshFunction


from anba4 import (
    material,
    anbax,
    DecoupleStiffness,
    PrincipalAxesRotationAngle,
    ComputeMassCenter,
    ComputeTensionCenter,
    ComputeShearCenter,
)
import argparse
import json
import yaml
import multiprocessing
import os
import numpy as np


def get_material_db(material_map):
    """
    Loads and processes a material database from a given material map file.

    Args:
        material_map (str): Path to the material map file. The file should be a JSON file
                            containing a mapping of material IDs to material properties.

    Returns:
        dict: A dictionary where keys are material IDs and values are material objects.
              The material objects can be either OrthotropicMaterial or IsotropicMaterial
              depending on the properties defined in the material database.

    Raises:
        AssertionError: If the material_map file does not exist.
        SystemExit: If the material map does not contain a link to a material database.

    Notes:
        - The material map file should contain a key "matdb" that points to the material
          database file (YAML format).
        - If the material database contains a material with ID "-1", it is assigned to
          the material map with the same ID.
        - The material properties for orthotropic materials are expected to include
          "tEx", "tEy", "tEz", "tGxy", "tGxz", "tGyz", "tnuxy", "tnuxz", and "tnuyz".
        - The material properties for isotropic materials are expected to include "E" or
          "Ex", and "nu".
        - The density of the material is expected to be defined by the "rho" key.
    """
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
    """
    Export unit strain fields to an XDMF file.

    This function computes and writes the strain fields for three unit directions
    ([1, 0, 0], [0, 1, 0], [0, 0, 1]) to an XDMF file. The strain fields are written
    at different time steps (0.0, 1.0, 2.0) to the specified result file.

    Parameters:
    anba (object): An object that contains the method `strain_field` and attribute `STRAIN`.
    result_file_name (str): The name of the XDMF file to write the results to.

    Returns:
    None
    """
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


def solve_anba4(mesh_filename, material_db_filename):
    """
    Solves the ANBA4 problem for a given mesh and material database.

    Parameters:
    mesh_filename (str): The filename of the mesh file in XDMF format.
    material_db_filename (str): The filename of the material database file.

    Returns:
    None

    The function performs the following steps:
    1. Checks if the mesh and material database files exist.
    2. Checks if the output JSON file already exists.
    3. Reads the mesh and material database.
    4. Initializes material parameters and orientations.
    5. Maps material indices from the mesh file.
    6. Sets material IDs and plane angles.
    7. Builds the material property library.
    8. Computes the stiffness and mass matrices.
    9. Decouples the stiffness matrix.
    10. Computes principal axes rotation angle, mass center, tension center, and shear center.
    11. Exports the results to an XDMF file.
    12. Writes the results to a JSON file.

    The output JSON file contains:
    - name: The name of the mesh file.
    - stiffness: The stiffness matrix.
    - mass_matrix: The mass matrix.
    - decoupled_stiffness: The decoupled stiffness matrix.
    - principal_axes_rotation: The principal axes rotation angle.
    - mass_center: The mass center.
    - tension_center: The tension center.
    - shear_center: The shear center.
    """
    if not os.path.isfile(mesh_filename):
        print(f"** Mesh file {mesh_filename} not found")
        return

    if not os.path.isfile(material_db_filename):
        print(f"** Material database file {material_db_filename} not found")
        return

    json_output_filename = f"{mesh_filename}.json"
    if os.path.isfile(json_output_filename):
        print(f"** Output file {json_output_filename} already exists, skipping")
        return

    print(f"run {mesh_filename}")
    material_db = get_material_db(material_db_filename)

    infile = XDMFFile(mesh_filename)

    mesh = Mesh()
    infile.read(mesh)

    pv_mesh = pv.read(mesh_filename)

    # Basic material parameters. 9 is needed for orthotropic materials.
    # TODO materials and orientations
    materials = MeshFunction("size_t", mesh, mesh.topology().dim())
    fiber_orientations = MeshFunction("double", mesh, mesh.topology().dim())
    plane_orientations = MeshFunction("double", mesh, mesh.topology().dim())

    # material indices (from material_map)
    material_unique_indices = np.unique(pv_mesh.cell_data["mat"])

    # map for materials in anba (index in material_unique_indices)
    material_map_0 = dict(
        zip(material_unique_indices, range(len(material_unique_indices)))
    )

    material_ids = [material_map_0[i] for i in pv_mesh.cell_data["mat"].tolist()]

    plane_angles = list(pv_mesh.cell_data["angle2"])

    materials.set_values(material_ids)

    # TODO, doesn't work for off axis laminates for now
    fiber_orientations.set_all(0.0)

    # transverse orientations
    plane_orientations.set_values(plane_angles)

    # Build material property library.

    material_library = [material_db[i] for i in material_unique_indices]

    anba = anbax(
        mesh, 2, material_library, materials, plane_orientations, fiber_orientations
    )
    stiffness = anba.compute()

    result_filename = mesh_filename.replace(".xdmf", "_results.xdmf")

    export_unit_strains(anba, result_filename)

    stiffness_matrix = stiffness.getDenseArray()

    mass = anba.inertia()

    mass_matrix = mass.getDenseArray()

    decoupled_stiffness = DecoupleStiffness(stiffness)

    principal_axes_rotation_angle = PrincipalAxesRotationAngle(decoupled_stiffness)

    mass_center = ComputeMassCenter(mass)
    tension_center = ComputeTensionCenter(stiffness_matrix)
    shear_center = ComputeShearCenter(stiffness_matrix)

    output = {
        "name": mesh_filename,
        "stiffness": stiffness_matrix.tolist(),
        "mass_matrix": mass_matrix.tolist(),
        "decoupled_stiffness": decoupled_stiffness.tolist(),
        "principal_axes_rotation": principal_axes_rotation_angle,
        "mass_center": mass_center,
        "tension_center": tension_center,
        "shear_center": shear_center,
    }

    with open(json_output_filename, "w") as write_file:
        json.dump(output, write_file, indent=4)


def main(args):
    """
    Main function to solve ANBA4 problems using multiprocessing.

    Args:
        args: An object containing the following attributes:
            - material_db: The material database to be used in the solver.
            - meshes: A list of mesh objects to be processed.

    The function utilizes all available CPU cores to parallelize the solving process.
    Each mesh is processed using the `solve_anba4` function with the provided material database.
    """
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.starmap(solve_anba4, [(mesh, args.material_db) for mesh in args.meshes])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="run a series of sections through anba4"
    )
    parser.add_argument("meshes", nargs="*", help="List of mesh files")
    parser.add_argument("material_db", help="Material map JSON file")
    # parser.add_argument("--debug", action="store_true", help="Run in debug mode")
    args = parser.parse_args()
    main(args)
