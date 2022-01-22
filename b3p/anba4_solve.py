#! /usr/bin/env python3
from dolfin import *
from anba4 import *
import argparse
import multiprocessing
import json


def run_mesh(meshname):
    mesh = Mesh()
    with XDMFFile(meshname) as infile:
        infile.read(mesh)

    # Basic material parameters. 9 is needed for orthotropic materials.
    # TODO materials and orientations
    E = 1.0
    nu = 0.33
    # Assmble into material mechanical property Matrix.

    matMechanicProp = [E, nu]
    # Meshing domain.
    # CompiledSubDomain
    materials = MeshFunction("size_t", mesh, mesh.topology().dim())
    fiber_orientations = MeshFunction("double", mesh, mesh.topology().dim())
    plane_orientations = MeshFunction("double", mesh, mesh.topology().dim())

    materials.set_all(0)
    fiber_orientations.set_all(0.0)
    plane_orientations.set_all(90.0)

    # Build material property library.
    mat1 = material.IsotropicMaterial(matMechanicProp, 1.0)

    matLibrary = []
    matLibrary.append(mat1)

    anba = anbax(mesh, 2, matLibrary, materials, plane_orientations, fiber_orientations)
    stiff = anba.compute()

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

    with open(meshname + ".json", "w") as write_file:
        json.dump(output, write_file, indent=4)


def main():
    p = argparse.ArgumentParser(description="run a series of sections through anba4")
    p.add_argument("meshes", nargs="*")

    args = p.parse_args()

    p = multiprocessing.Pool()
    r = p.map_async(
        run_mesh, args.meshes
    )  # run async seems to avoid a sporadic blocking error
    r.wait()


if __name__ == "__main__":
    main()
