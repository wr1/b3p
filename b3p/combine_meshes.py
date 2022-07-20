#! /usr/bin/env python3

# import vtk
import argparse
import numpy as np
import pyvista as pv
import time
import multiprocessing

# meshes = []


# def add_points(inp):
#     global meshes
#     # print(len(meshes), inp[1])
#     if inp[0] == 0:
#         # print(inp[:2])
#         meshes[inp[1]].point_data[inp[2]] = inp[3]
#     elif inp[0] == 1:
#         print(inp[:3])
#         meshes[inp[1]].cell_data[inp[2]] = inp[3]


def main():
    # global meshes
    p = argparse.ArgumentParser(
        description="Join a series of meshes, i.e. a shell and n web meshes together into a single vtu"
    )
    p.add_argument("meshes", nargs="*")
    p.add_argument("--out", default="__joined_mesh.vtu", help="output file name")
    args = p.parse_args()

    meshes = []
    for i in args.meshes:
        meshes.append(pv.read(i))

    all_pd = [
        (j, x.point_data[j].shape, x.point_data[j].dtype)
        for x in meshes
        for j in x.point_data.keys()
    ]
    all_cd = [
        (j, x.cell_data[j].shape, x.cell_data[j].dtype)
        for x in meshes
        for j in x.cell_data.keys()
    ]

    tic = time.time()

    join = []

    for i in range(len(meshes)):
        m = meshes[i]
        for j in all_pd:
            if j[0] not in m.point_data:
                a = np.zeros((m.n_points, j[1][1] if len(j[1]) > 1 else 1), dtype=j[2])
                m.point_data[j[0]] = a
                # join.append((0, i, j[0], a))
        for j in all_cd:
            if j[0] not in m.cell_data:
                a = np.zeros((m.n_cells, j[1][1] if len(j[1]) > 1 else 1), dtype=j[2])
                m.cell_data[j[0]] = a
                # join.append((1, i, j[0], a))

    # pool = multiprocessing.Pool()
    # pool.map(add_points, join)
    # pool.join()

    toc = time.time()

    out = meshes[0].merge(meshes[1:])

    toc2 = time.time()

    print(toc - tic, toc2 - toc)

    out.save(args.out)
    print("written mesh to %s" % args.out)

    """
    append = vtk.vtkAppendFilter()
    append.MergePointsOn()
    append.ToleranceIsAbsoluteOn()
    append.SetTolerance(1e-3)

    allcellarrays = []
    meshes = []
    # loop over all meshes and get associated arrays
    for i in args.meshes:
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(i)
        reader.Update()
        meshes.append(reader.GetOutput())

        for j in range(meshes[-1].GetCellData().GetNumberOfArrays()):
            allcellarrays.append(
                (
                    meshes[-1].GetCellData().GetArrayName(j),
                    meshes[-1].GetCellData().GetArray(j).GetNumberOfComponents(),
                )
            )

    # make sure all meshes have all arrays (add zero arrays) so that they show up after the merge
    for i in meshes:
        for j in allcellarrays:
            if not i.GetCellData().HasArray(j[0]):
                arr = vtk.vtkFloatArray()
                arr.SetNumberOfComponents(j[1])
                arr.SetNumberOfTuples(i.GetNumberOfCells())
                for k in range(j[1]):
                    arr.FillComponent(k, 0.0)
                arr.SetName(j[0])
                i.GetCellData().AddArray(arr)

    for i in meshes:
        append.AddInputData(i)
    append.Update()
    append = append.GetOutput()

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(args.out)
    writer.SetInputData(append)
    writer.Update()
    writer.Write()
    print("written mesh to %s" % args.out)"""


if __name__ == "__main__":
    main()
