#! /usr/bin/env python3

import pyvista as pv
import numpy as np
import multiprocessing
import vtk
import argparse
import time
import json
import os
import yaml


def make_shell_section(elem_id, plyarray, merge_adjacent_plies=True, zero_angle=True):
    plies = []
    # if material in next ply is the same, add it to the previous ply,
    # NOTE this does not take ply angle into account so it only works
    # for multidirectional fibre mats at angle==0 for now
    for j in plyarray[np.where(plyarray[:, 1] > 1e-6)]:
        if plies and plies[-1][1] == j[0] and merge_adjacent_plies:
            plies[-1][0] += j[1] * 1e-3
        else:
            plies.append([j[1] * 1e-3, j[0]])
    if zero_angle:
        return "".join("%f,,m%i,0\n" % tuple(i) for i in plies)
    else:
        return "".join("%f,,m%i,or%i\n" % tuple(i + [elem_id + 1]) for i in plies)


def material_db_to_ccx(materials, matmap=None, force_iso=False):
    "find the material db and write properties to a ccx block"
    # 2 files are relevant, a material map that maps the material ID (integer)
    # in the VTK file to a key (string) in the material database and the material
    # database itself
    mat_db = None
    if os.path.isfile(matmap):  # check if the material map file is there
        gdir = os.path.dirname(matmap)
        mm = json.load(open(matmap, "r"))
        if "matdb" in mm:  # check if the material map file points to a material db
            mat_db = yaml.load(
                open(os.path.join(gdir, mm["matdb"])), Loader=yaml.CLoader
            )
        else:
            exit(
                "material map available, but no link to material db, need matdb definition to do FEA"
            )
    else:
        exit("no material map defined")

    mm_inv = {v: k for k, v in mm.items()}

    matblock = ""
    for i in materials:
        if i > 1e-6:
            material_properties = mat_db[mm_inv[int(i)]]

            # print(material_properties)
            matblock += f'** material: {material_properties["name"]}\n'
            # matblock += f"** {str(material_properties)}\n"

            if (
                "vf" in material_properties
                or "C" in material_properties
                and not force_iso
            ):
                print(material_properties["name"], "is assumed to be orthotropic")
                C = np.array(material_properties["C"])
                # https://github.com/rsmith-nl/lamprop/blob/410ebfef2e14d7cc2988489ca2c31103056da38f/lp/text.py#L96
                # https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node193.html
                matblock += "** orthotropic material\n"
                matblock += "*material,name=m%i\n*elastic,type=ortho\n" % i
                D = C
                D[0, 3] = C[0, 5]
                D[0, 5] = C[0, 3]
                D[1, 3] = C[1, 5]
                D[1, 5] = C[1, 3]
                D[2, 3] = C[2, 5]
                D[2, 5] = C[2, 3]
                D[3, 3] = C[5, 5]
                D[5, 5] = C[3, 3]
                D *= 1e6
                matblock += (
                    f"{D[0,0]:.4g},{D[0,1]:.4g},{D[1,1]:.4g},"
                    + f"{D[0,2]:.4g},{D[1,2]:.4g},{D[2,2]:.4g},"
                    + f"{D[3,3]:.4g},{D[4,4]:.4g},\n"
                    + f"{D[5,5]:.4g},293\n"
                )
            elif "e11" in material_properties and not force_iso:
                print(material_properties["name"], "has engineering constants")

                matblock += "** orthotropic material\n"
                matblock += (
                    "*material,name=m%i\n*elastic,type=engineering constants\n" % i
                )
                matblock += (
                    f"{material_properties['e11']:.4g},{material_properties['e22']:.4g},{material_properties['e33']:.4g},"
                    + f"{material_properties['nu12']:.4g},{material_properties['nu31']:.4g},{material_properties['nu23']:.4g},"
                    + f"{material_properties['g12']:.4g},{material_properties['g31']:.4g},\n"
                    + f"{material_properties['g23']:.4g},293\n"
                )

            else:
                print(material_properties["name"], "is assumed to be isotropic")
                print(material_properties)
                nu = min(
                    0.45,
                    max(
                        0.1,
                        float(material_properties["nu"])
                        if "nu" in material_properties
                        else material_properties["nu12"],
                    ),
                )  # min()
                E = 1e6 * float(
                    material_properties["tEx"]
                    if "tEx" in material_properties
                    else material_properties["E"]
                )
                matblock += "** isotropic material\n"
                matblock += "*material,name=m%i\n*elastic,type=iso\n" % i
                matblock += f"{E:.4g},{nu:.4g},293\n"

    return matblock


def format_eset(name, eids, prefix):
    out = f"*elset,elset={name.replace(prefix,'').lstrip('0')}\n"
    for i in range(len(eids)):
        out += f"{eids[i]}"
        out += "\n" if (i % 16 == 15) else ","
    if out[-1] == ",":
        out = out[:-1] + "\n"
    return out


def compute_groups(grid, prefix):
    gr = ""
    for i in grid.cell_data:
        if i.startswith(prefix):
            if len(grid.cell_data[i].shape) == 1:
                eids = np.where(grid.cell_data[i] > 0)[0] + 1
            elif len(grid.cell_data[i].shape):
                eids = np.where(grid.cell_data[i][:, 1] > 0)[0] + 1

            gr += format_eset(i, eids, prefix)

    return gr


def nodebuffer(grid):
    buf = "*node,nset=nall\n"
    for n, i in enumerate(grid.points):
        buf += "%i,%f,%f,%f\n" % tuple([n + 1] + list(i))
    return buf


def element_buffer_quadratic(grid):
    conn = grid.cell_connectivity.reshape(
        (
            grid.GetNumberOfCells(),
            int(grid.cell_connectivity.shape[0] / grid.GetNumberOfCells()),
        )
    )
    buf = ""
    # export the elements
    for n, i in enumerate(conn):
        buf += "*element,type=s8r,elset=e%i\n" % (n + 1)
        buf += "%i,%i,%i,%i,%i,%i,%i,%i,%i\n" % tuple([n + 1] + list(i + 1))

    return buf


def element_buffer_linear(grid):
    conn = grid.cell_connectivity.reshape(
        (
            grid.GetNumberOfCells(),
            int(grid.cell_connectivity.shape[0] / grid.GetNumberOfCells()),
        )
    )
    buf = ""
    # export the elements
    for n, i in enumerate(conn):
        buf += "*element,type=s4,elset=e%i\n" % (n + 1)
        buf += "%i,%i,%i,%i,%i\n" % tuple([n + 1] + list(i + 1))

    return buf


def orientation_buffer(grid, add_centers=False):
    buf = ""
    # write orientation TODO match with element orientation, for now just align with z-axis
    for n in range(grid.GetNumberOfCells()):
        xdir, ydir = grid.cell_data["x_dir"][n], grid.cell_data["y_dir"][n]
        center = grid.cell_data["centers"][n]
        buf += "*orientation,name=or%i,system=rectangular\n" % (n + 1)
        buf += "%.4g,%.4g,%.4g,%.4g,%.4g,%.4g" % tuple(xdir.tolist() + ydir.tolist())
        buf += "," + ",".join(center.astype(str)) + "\n" if add_centers else "\n"
    return buf


def get_loadcases(mesh, multiplier=1.0):
    loadcases = {}
    for i in mesh.point_data:
        if i.startswith("lc_"):
            print(f"loadcase {i}")
            # forces are interpolated to midside nodes, causing the sum of forces to be off,
            # compute a multiplier from the sum of the forces in the linear model here
            multiplier = 1.0  # TODO fix for quadratic meshes # mesh.point_data[i].sum() / mesh.point_data[i].sum()
            lbuf = f"** {i}\n*step\n*static\n*cload\n"
            ld = mesh.point_data[i] * multiplier
            for n, j in enumerate(ld):
                if j[0] ** 2 > 1e-8:
                    lbuf += "%i,1,%f\n" % (n + 1, j[0])
                if j[1] ** 2 > 1e-8:
                    lbuf += "%i,2,%f\n" % (n + 1, j[1])

            root = np.where(mesh.points[:, 2] == mesh.points[:, 2].min())
            lbuf += "*boundary,op=new\n"
            for j in root[0]:
                lbuf += "%i,1,3\n" % (j + 1)

            lbuf += "*node file,output=3d\nU,RF\n*EL FILE\nS,E\n*node print,nset=nall\nrf\n*end step\n"

            loadcases[i] = lbuf

    return loadcases


def main():
    p = argparse.ArgumentParser(description="Translate a blade vtk mesh into ccx")
    p.add_argument(
        "grid",
        help="Grid file (vtu) including shear web(s) and shell (assuming you want to simulate a blade)",
    )
    p.add_argument("--out", default="test.inp", help="Output file name (input for ccx)")
    p.add_argument(
        "--matmap",
        default="temp/material_map.json",
        help="Mapping between material name and id in the vtu file",
    )
    p.add_argument(
        "--merge_adjacent_layers",
        action="store_true",
        default=False,
        help="Merge adjacent layers with same material",
    )
    p.add_argument(
        "--zeroangle",
        action="store_true",
        default=False,
        help="Zero angle set for plies, rather than reference to coordinate system. This can set for abaqus and hyperworks",
    )
    p.add_argument(
        "--single",
        action="store_true",
        default=False,
        help="If set, export all loadcases as *STEPs in a single file rather than making a separate file for each loadcase",
    )
    p.add_argument(
        "--quadratic", action="store_true", default=False, help="Use quadratic elements"
    )
    p.add_argument(
        "--add_centers",
        action="store_true",
        default=False,
        help="Add centers to orientation",
    )
    p.add_argument("--force_iso", action="store_true", default=False, help="Force iso")
    args = p.parse_args()

    grid = args.grid

    # order = 1

    print(pv.read(grid).points.shape)

    gr = pv.read(grid).threshold(value=(1e-6, 1e9), scalars="thickness")

    print(gr.points.shape)
    gr.cell_data["centers"] = gr.cell_centers().points

    if args.quadratic:
        lf = vtk.vtkLinearToQuadraticCellsFilter()
        lf.SetInputData(gr)
        lf.Update()
        quad = lf.GetOutput()
        mesh = pv.UnstructuredGrid(quad)
        # mesh = mesh
    else:
        mesh = pv.UnstructuredGrid(gr)

    # export the nodes
    buf = "*node,nset=nall\n"
    for n, i in enumerate(mesh.points):
        buf += "%i,%f,%f,%f\n" % tuple([n + 1] + list(i))

    buf = nodebuffer(mesh)

    if args.quadratic:  # order == 1:
        buf += element_buffer_quadratic(mesh)
    else:
        buf += element_buffer_linear(mesh)

    # else:
    #     raise ValueError("order must be 1 or 2")
    print(mesh.GetNumberOfCells())

    buf += "*elset,elset=Eall,GENERATE\n%i,%i\n" % (1, mesh.GetNumberOfCells())

    buf += compute_groups(mesh, "slab_thickness_")

    buf += compute_groups(mesh, "ply_")

    plykeys = [i for i in mesh.cell_data if i.startswith("ply_")]

    plydat = np.stack(mesh.cell_data[i] for i in plykeys)
    # get all materials of all plies
    materials = np.unique(plydat[:, :, 0])

    buf += orientation_buffer(mesh, args.add_centers)

    matblock = material_db_to_ccx(
        materials, matmap=args.matmap, force_iso=args.force_iso
    )

    # exit()
    buf += matblock

    tic = time.perf_counter()

    # p = multiprocessing.Pool()

    # blx = p.map(
    #     make_shell_section,
    #     [
    #         (i, plydat[:, i, :], args.merge_adjacent_layers)
    #         for i in range(plydat.shape[1])
    #     ],
    # )

    blx = [
        make_shell_section(
            i, plydat[:, i, :], args.merge_adjacent_layers, args.zeroangle
        )
        for i in range(plydat.shape[1])
    ]

    toc = time.perf_counter()
    print("time spent creating shell sections %f" % (toc - tic))

    comps = "".join(
        "*shell section,composite,elset=e%i,offset=-.5\n%s" % (n + 1, i)
        for n, i in enumerate(blx)
    )
    buf += comps

    loadcases = get_loadcases(mesh)

    # print(loadcases)

    # write a full ccx file for each loadcase, assuming parallel execution
    if args.single:
        output = buf + "".join(loadcases.values())
        open(args.out, "w").write(output)
        print(f"written ccx input file with all loadcases to {args.out}")
    else:
        for i in loadcases:
            output = buf + loadcases[i]
            of = args.out.replace(".inp", f"_{i}.inp")
            open(of, "w").write(output)
            print(f"written ccx input file to {of}")


if __name__ == "__main__":
    main()
