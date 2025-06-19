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

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def get_material_db(material_map, unit_factor=1):
    """Load and process a material database from a JSON material map file."""
    if not os.path.isfile(material_map):
        raise FileNotFoundError(f"Material map file not found: {material_map}")

    with open(material_map, "r") as f:
        mm1 = json.load(f)
        mm = mm1.get("map", mm1)  # Support both old and new format
        mat_db = mm1.get("matdb", None)

    if "-1" in mat_db:
        mm["-1"] = -1

    materials = {}
    for mat_name, mesh_id in mm.items():
        if mat_name == "matdb":
            continue
        matdb_id = mat_name
        if matdb_id not in mat_db:
            logger.error(
                f"Material ID {matdb_id} from material_map not found in materials.yaml"
            )
            raise KeyError(f"Material ID {matdb_id} not found in materials.yaml")
        matdb_entry = mat_db[matdb_id]
        logger.debug(
            f"Processing material {matdb_id} (mesh ID {mesh_id}): {matdb_entry}"
        )
        density = matdb_entry.get("rho", matdb_entry.get("density", 1.0))

        if "e11" in matdb_entry:
            required_keys = [
                "e11",
                "e22",
                "e33",
                "g12",
                "g13",
                "g23",
                "nu12",
                "nu13",
                "nu23",
            ]
            if not all(k in matdb_entry for k in required_keys):
                raise ValueError(
                    f"Missing orthotropic properties for material {matdb_id}: {required_keys}"
                )
            matMechanicProp = np.zeros((3, 3))
            matMechanicProp[0, 2] = matdb_entry["e11"] * unit_factor
            matMechanicProp[0, 1] = matdb_entry["e22"] * unit_factor
            matMechanicProp[0, 0] = matdb_entry["e33"] * unit_factor
            matMechanicProp[1, 2] = matdb_entry["g23"] * unit_factor
            matMechanicProp[1, 1] = matdb_entry["g13"] * unit_factor
            matMechanicProp[1, 0] = matdb_entry["g12"] * unit_factor
            matMechanicProp[2, 2] = matdb_entry["nu23"]
            matMechanicProp[2, 1] = matdb_entry["nu13"]
            matMechanicProp[2, 0] = matdb_entry["nu12"]
            materials[matdb_id] = material.OrthotropicMaterial(matMechanicProp, density)

        else:
            if "E" not in matdb_entry and "Ex" not in matdb_entry:
                raise ValueError(
                    f"Missing 'E' or 'Ex' for isotropic material {matdb_id}. Entry: {matdb_entry}"
                )
            if "nu" not in matdb_entry:
                raise ValueError(
                    f"Missing 'nu' for isotropic material {matdb_id}. Entry: {matdb_entry}"
                )
            E = matdb_entry.get("E", matdb_entry.get("Ex")) * unit_factor
            nu = min(matdb_entry["nu"], 0.49)
            materials[matdb_id] = material.IsotropicMaterial([E, nu], density)

    mm_inv = {str(v): k for k, v in mm.items() if k != "matdb"}
    logger.debug(f"Material map inverse: {mm_inv}")

    return materials, mm_inv


def export_unit_strains(anba, result_file_name, pv_mesh):
    """Export strain fields for unit displacements to separate VTK files."""

    for disp, suffix in [([1, 0, 0], "mx"), ([0, 1, 0], "my"), ([0, 0, 1], "mz")]:
        anba.strain_field([0, 0, 0], disp, "local", "paraview")
        nodevalues = anba.STRAIN.vector().get_local().reshape(-1, 6)
        pv_mesh.cell_data[f"strain_{suffix}"] = nodevalues

    pv_mesh.save(result_file_name)


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
    material_db, mm_inv = get_material_db(material_db_filename)

    with XDMFFile(mesh_filename) as infile:
        mesh = Mesh()
        infile.read(mesh)
    # Read .xdmf directly with PyVista, handling potential lack of time steps
    pv_mesh = pv.read(mesh_filename.replace(".xdmf", "_tri.vtp"))

    materials = MeshFunction("size_t", mesh, mesh.topology().dim())
    fiber_orientations = MeshFunction("double", mesh, mesh.topology().dim())
    plane_orientations = MeshFunction("double", mesh, mesh.topology().dim())

    material_unique_indices = np.unique(pv_mesh.cell_data["mat"])
    material_unique_indices.sort()
    logger.debug(f"Mesh material IDs: {material_unique_indices.tolist()}")

    material_map = {}
    for mesh_id in material_unique_indices:
        mesh_id_str = str(mesh_id)
        if mesh_id_str not in mm_inv:
            logger.error(
                f"Mesh material ID {mesh_id_str} not found in material map inverse: {mm_inv}"
            )
            raise KeyError(f"Mesh material ID {mesh_id_str} not found in material map")
        material_map[mesh_id] = mm_inv[mesh_id_str]
    logger.debug(f"Material mapping: {material_map}")

    anba_map = {mesh_id: idx for idx, mesh_id in enumerate(material_unique_indices)}
    logger.debug(f"ANBA4 mapping (mesh ID to index): {anba_map}")

    material_ids = [anba_map[i] for i in pv_mesh.cell_data["mat"].tolist()]
    logger.debug(f"ANBA4 material IDs: {np.unique(material_ids).tolist()}")

    material_library = [
        material_db[material_map[mesh_id]] for mesh_id in material_unique_indices
    ]
    logger.debug(
        f"Material library order: {[material_map[mesh_id] for mesh_id in material_unique_indices]}"
    )

    materials.set_values(material_ids)
    fiber_orientations.set_all(0.0)  # TODO: Support off-axis laminates
    plane_orientations.set_values(pv_mesh.cell_data["angle2"].tolist())

    anba = anbax(
        mesh, 1, material_library, materials, plane_orientations, fiber_orientations
    )
    stiffness = anba.compute()
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
    with open(json_output_filename, "w") as f:
        json.dump(output, f, indent=4)

    result_file = mesh_filename.replace(".xdmf", "_results.vtp")
    export_unit_strains(anba, result_file, pv_mesh)
    # import pyvista as pv
    # mesh = pv.read(result_filename)
    # mesh.save(result_filename.replace(".xdmf", ".vtu"))
    logger.info(f"Wrote results to {json_output_filename} and {result_file}")


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
