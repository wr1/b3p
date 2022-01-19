#! /usr/bin/env python3

import vtk
import argparse
import numpy as np


def main():
    p = argparse.ArgumentParser(
        description="Join a series of meshes, i.e. a shell and n web meshes together into a single vtu"
    )
    p.add_argument("meshes", nargs="*")
    p.add_argument("--out", default="__joined_mesh.vtu", help="output file name ")
    args = p.parse_args()

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
    print("written mesh to %s" % args.out)


if __name__ == "__main__":
    main()