#! /usr/bin/env python3
from dolfin import *
from anba4 import *
import argparse
import multiprocessing
import json
import numpy as np
import yaml
import os 
import pyvista as pv 

def get_material_db(material_map):
    assert os.path.isfile(material_map)
    mm = json.load(open(material_map,'r'))
    gdir = os.path.dirname(material_map)
    
    mat_db = None
    if "matdb" in mm:  # check if the material map file points to a material db
        mat_db = yaml.load(
            open(os.path.join(gdir, mm["matdb"])), Loader=yaml.CLoader
        )
    else:
        exit(
            "material map available, but no link to material db, need matdb definition to do FEA"
        )


    mm_inv = {v: k for k, v in mm.items()}

    for i in mm_inv:
        if i !=mm['matdb']:      
            material = mat_db[mm_inv[i]]
            if 'vf' in material: # ortho material
                matMechanicProp = np.zeros((3,3))
                matMechanicProp[0,0] = material['tEx']*1e6 # e_xx
                matMechanicProp[0,1] = material['tEy']*1e6 # e_yy
                matMechanicProp[0,2] = material['tEz']*1e6 # e_zz
                matMechanicProp[1,0] = material['tGyz']*1e6#g_yz
                matMechanicProp[1,1] = material['tGxz']*1e6#g_xz
                matMechanicProp[1,2] = material['tGxy']*1e6#g_xy
                matMechanicProp[2,0] = material['tnuyz']#nu_zy
                matMechanicProp[2,1] = material['tnuxz']#nu_zx
                matMechanicProp[2,2] = material['tnuxy']#nu_xy
            else:
                
            print(i)






    return ( mm_inv, mat_db )





def run_mesh(meshname):

    infile = XDMFFile(meshname)
    mesh = Mesh()
    infile.read(mesh)

    pvmesh = pv.read(meshname)


    matid = MeshFunction("size_t", mesh, mesh.topology().dim())



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

    print(dir(materials))

    matuniq = np.unique(pvmesh.cell_data['mat'])


    # map for materials in this section 
    mat_map_0 = dict(zip(matuniq,range(len(matuniq)))) 

    materials.set_values([mat_map_0[i] for i in pvmesh.cell_data['mat'].tolist()])
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
    p.add_argument("matdb",help="material map json file")
    p.add_argument('--debug',action='store_true')
    args = p.parse_args()

    mdb = get_material_db(args.matdb)
    exit()
    if args.debug:
        for i in args.meshes:
            run_mesh(i)
    else:
        p = multiprocessing.Pool()
        r = p.map_async(
            run_mesh, args.meshes
        )  # run async seems to avoid a sporadic blocking error
        r.wait()


if __name__ == "__main__":
    main()
