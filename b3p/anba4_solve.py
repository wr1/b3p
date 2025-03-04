#!/usr/bin/env python3
import argparse
import json
import logging
import multiprocessing
import os
import traceback
import numpy as np
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
import yaml

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def get_material_db(material_map, unit_factor=1):
    """
    Load and process a material database from a JSON material map file.

    Args:
        material_map (str): Path to the JSON material map file.
        unit_factor (float): Conversion factor for material properties (default: 1e6, e.g., MPa to Pa).

    Returns:
        dict: Mapping of material IDs to material objects (OrthotropicMaterial or IsotropicMaterial).

    Raises:
        FileNotFoundError: If the material map file or material database file doesn't exist.
        ValueError: If required material properties are missing or invalid.
    """
    if not os.path.isfile(material_map):
        raise FileNotFoundError(f"Material map file not found: {material_map}")

    with open(material_map, "r") as f:
        mm = json.load(f)

    gdir = os.path.dirname(material_map)
    if "matdb" not in mm:
        raise ValueError(
            "Material map must contain 'matdb' key pointing to material database"
        )

    mat_db_path = os.path.join(gdir, mm["matdb"])
    if not os.path.isfile(mat_db_path):
        raise FileNotFoundError(f"Material database file not found: {mat_db_path}")

    with open(mat_db_path, "r") as f:
        mat_db = yaml.safe_load(f)

    if "-1" in mat_db:
        mm["-1"] = -1

    mm_inv = {v: k for k, v in mm.items() if v != mm["matdb"]}
    materials = {}
    shadow_materials = {}

    for i in mm_inv:
        matdb_entry = mat_db[mm_inv[i]]
        logger.debug(f"Processing material {mm_inv[i]}: {matdb_entry}")
        density = matdb_entry.get("rho", matdb_entry.get("density", 1.0))

        if "e11" in matdb_entry:  # Orthotropic material
            # required_keys = [
            #     "e11",
            #     "e22",
            #     "e33",
            #     "g12",
            #     "g13",
            #     "g23",
            #     "nu12",
            #     "nu13",
            #     "nu23",
            # ]
            # if not all(k in matdb_entry for k in required_keys):
            #     raise ValueError(
            #         f"Missing orthotropic properties for material {mm_inv[i]}: {required_keys}"
            #     )

            # matMechanicProp = np.zeros((3, 3))
            # matMechanicProp[0, 0] = matdb_entry["e11"] * unit_factor  # E_xx
            # matMechanicProp[0, 1] = matdb_entry["e22"] * unit_factor  # E_yy
            # matMechanicProp[0, 2] = matdb_entry["e33"] * unit_factor  # E_zz
            # matMechanicProp[1, 0] = matdb_entry["g23"] * unit_factor  # G_xy
            # matMechanicProp[1, 1] = matdb_entry["g13"] * unit_factor  # G_xz
            # matMechanicProp[1, 2] = matdb_entry["g12"] * unit_factor  # G_yz
            # matMechanicProp[2, 0] = matdb_entry["nu23"]  # nu_xy
            # matMechanicProp[2, 1] = matdb_entry["nu13"]  # nu_xz
            # matMechanicProp[2, 2] = matdb_entry["nu12"]  # nu_yz
            # materials[i] = material.OrthotropicMaterial(matMechanicProp, density)

            # Simplified orthotropic material
            E = matdb_entry["e11"] * unit_factor
            nu = matdb_entry["nu12"]
            materials[i] = material.IsotropicMaterial([E, nu], density)
            shadow_materials[i] = [E, nu]
        else:  # Isotropic material
            if "E" not in matdb_entry and "Ex" not in matdb_entry:
                raise ValueError(
                    f"Missing 'E' or 'Ex' for isotropic material {mm_inv[i]}. Entry: {matdb_entry}"
                )
            if "nu" not in matdb_entry:
                raise ValueError(
                    f"Missing 'nu' for isotropic material {mm_inv[i]}. Entry: {matdb_entry}"
                )
            E = matdb_entry.get("E", matdb_entry.get("Ex")) * unit_factor
            nu = min(matdb_entry["nu"], 0.49)
            materials[i] = material.IsotropicMaterial([E, nu], density)
            shadow_materials[i] = [E, nu]
    return materials, shadow_materials


def export_unit_strains(anba, result_file_name):
    """Export strain fields for unit displacements to an XDMF file."""
    result_file = XDMFFile(result_file_name)
    result_file.parameters["functions_share_mesh"] = True
    result_file.parameters["rewrite_function_mesh"] = False
    result_file.parameters["flush_output"] = True

    anba.strain_field([0, 0, 0], [1, 0, 0], "local", "paraview")
    result_file.write(anba.STRAIN, t=0.0)
    anba.strain_field([0, 0, 0], [0, 1, 0], "local", "paraview")
    result_file.write(anba.STRAIN, t=1.0)
    anba.strain_field([0, 0, 0], [0, 0, 1], "local", "paraview")
    result_file.write(anba.STRAIN, t=2.0)

    logger.info(f"Wrote strain results to {result_file_name}")


def solve_anba4(mesh_filename, material_db_filename):
    """Solve the ANBA4 problem for a mesh and material database."""
    if not mesh_filename.endswith(".xdmf"):
        logger.error(f"Expected .xdmf mesh file, got {mesh_filename}")
        return
    if not os.path.isfile(mesh_filename):
        logger.error(f"Mesh file not found: {mesh_filename}")
        return
    if not os.path.isfile(material_db_filename):
        logger.error(f"Material database file not found: {material_db_filename}")
        return

    json_output_filename = f"{mesh_filename}.json"
    if os.path.isfile(json_output_filename):
        logger.info(f"Output file {json_output_filename} exists, skipping")
        return

    logger.info(f"Processing {mesh_filename}")
    material_db, shadow = get_material_db(material_db_filename)

    with XDMFFile(mesh_filename) as infile:
        mesh = Mesh()
        infile.read(mesh)
    pv_mesh = pv.read(mesh_filename.replace(".xdmf", ".vtp"))

    materials = MeshFunction("size_t", mesh, mesh.topology().dim())
    fiber_orientations = MeshFunction("double", mesh, mesh.topology().dim())
    plane_orientations = MeshFunction("double", mesh, mesh.topology().dim())

    # Get unique material indices from the mesh
    material_unique_indices = np.unique(pv_mesh.cell_data["mat"])
    material_unique_indices.sort()  # Sort if order matters

    # Create a mapping from original material IDs to a contiguous range
    material_map = dict(
        zip(material_unique_indices, range(len(material_unique_indices)))
    )
    logger.debug(f"Material mapping: {material_map}")

    # Remap material IDs from the mesh to the contiguous range
    material_ids = [material_map[i] for i in pv_mesh.cell_data["mat"].tolist()]
    plane_angles = pv_mesh.cell_data["angle2"].tolist()

    # Set the remapped material IDs into the mesh function
    materials.set_values(material_ids)
    fiber_orientations.set_all(0.0)  # TODO: Support off-axis laminates
    plane_orientations.set_values(plane_angles)

    # Build material library in the order of remapped indices
    material_library = []
    for original_id in material_unique_indices:
        if original_id not in material_db:
            logger.error(f"Material ID {original_id} not found in material database")
            raise KeyError(f"Material ID {original_id} not found in material database")
        material_library.append(material_db[original_id])

    # Verify alignment
    if len(material_library) != len(material_unique_indices):
        logger.error(
            f"Material library size ({len(material_library)}) does not match unique indices ({len(material_unique_indices)})"
        )
        raise ValueError("Mismatch between material library and unique indices")

    print(material_ids[392], shadow[material_ids[392]])

    # Debug: Check a specific material (e.g., index 392)
    # if len(material_ids) > 392:
    #     remapped_id = material_ids[392]
    #     logger.debug(
    #         f"Material ID at index 392: {remapped_id}, Material properties: {dir(material_library[remapped_id])}"
    #     )

    # Print the E modulus of the material at index 392
    # if len(material_ids) > 392:
    #     remapped_id = material_ids[392]
    #     material_obj = material_library[remapped_id]
    #     if isinstance(material_obj, material.IsotropicMaterial):
    #         # For isotropic materials, E is the first element of the property list
    #         e_modulus = material_obj.mechanic_prop[
    #             0
    #         ]  # Assuming mechanic_prop holds [E, nu]
    #         logger.info(f"Material at index 392 (ID {remapped_id}) has E = {e_modulus}")
    #     elif isinstance(material_obj, material.OrthotropicMaterial):
    #         # For orthotropic materials, e11 is the first modulus (xx direction)
    #         e_modulus = material_obj.mechanic_prop[0, 0]
    #         # Assuming mechanic_prop is a 3x3 array
    #         logger.info(
    #             f"Material at index 392 (ID {remapped_id}) has e11 = {e_modulus}"
    #         )
    #     else:
    #         logger.warning(f"Unknown material type at index 392: {type(material_obj)}")

    # Initialize ANBA4 solver
    anba = anbax(
        mesh, 1, material_library, materials, plane_orientations, fiber_orientations
    )
    stiffness = anba.compute()
    stiffness_matrix = stiffness.getDenseArray()

    mass = anba.inertia()
    mass_matrix = mass.getDenseArray()

    # Compute additional properties
    decoupled_stiffness = DecoupleStiffness(stiffness)
    principal_axes_rotation_angle = PrincipalAxesRotationAngle(decoupled_stiffness)
    mass_center = ComputeMassCenter(mass)
    tension_center = ComputeTensionCenter(stiffness_matrix)
    shear_center = ComputeShearCenter(stiffness_matrix)

    # Save results
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
    with open(json_output_filename, "w") as f:
        json.dump(output, f, indent=4)

    result_filename = mesh_filename.replace(".xdmf", "_results.xdmf")
    export_unit_strains(anba, result_filename)


def run_solver(mesh, material_db):
    """Wrapper to run solve_anba4 with error handling for multiprocessing."""
    try:
        solve_anba4(mesh, material_db)
    except Exception as e:
        logger.error(f"Failed to process {mesh}: {str(e)}\n{traceback.format_exc()}")


def main(args):
    """Process multiple meshes using ANBA4 solver in parallel."""
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.starmap(run_solver, [(mesh, args.material_db) for mesh in args.meshes])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run ANBA4 solver on a series of mesh files."
    )
    parser.add_argument("meshes", nargs="*", help="List of XDMF mesh files")
    parser.add_argument("material_db", help="Material map JSON file")
    args = parser.parse_args()
    main(args)
