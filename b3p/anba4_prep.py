#! /usr/bin/env python3
import pyvista as pv
import argparse
import multiprocessing


def conv3d_2d(mesh):
    """hack to translate a mesh to 2D format, which is what fenics/anba needs
    :param mesh: path to the mesh file
    """
    x = open(mesh, "r").read()
    nb = x[: x.find("Topology")]
    rest = x[x.find("Topology") :]
    nb = nb.replace("0.0000000e+00\n", "\n")
    nb = nb.replace(' 3"', ' 2"')
    nb = nb.replace("XYZ", "XY")
    open(mesh, "w").write(nb + rest)


def vtp2xdmf(vtp):
    """convert a vtp file to xdmf format and translate to 2D"""
    assert vtp.endswith(".vtp")
    mesh = pv.read(vtp)
    tri = mesh.triangulate()  # fenics doesn't do mixed element types
    tri.points[:, 2] = 0
    xd = vtp.replace(".vtp", ".xdmf")
    pv.save_meshio(xd, tri, data_format="XML")
    conv3d_2d(xd)
    print(f"converted {vtp} to {xd}")


def anba4_prep(section_meshes):
    p = multiprocessing.Pool()
    p.map(vtp2xdmf, section_meshes)
    p.close()


def main():
    p = argparse.ArgumentParser(
        description="translate section meshes from vtk to XDMF and 2D"
    )
    p.add_argument("sections", nargs="*", help="section meshes in .vtp format")
    args = p.parse_args()
    anba4_prep(args.sections)


if __name__ == "__main__":
    main()
