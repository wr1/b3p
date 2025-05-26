#! /usr/bin/env python3
import pyvista as pv
import argparse
import multiprocessing
import os
import logging

logger = logging.getLogger(__name__)


def conv3d_2d(mesh):
    """Hack to translate a 3D mesh to 2D format, which is required by fenics/anba.

    This function reads a mesh file, modifies its content to convert it from
    3D to 2D format, and writes the modified content back to the file.

    Parameters:
    mesh (str): The path to the mesh file.

    The function performs the following modifications:
    - Removes occurrences of "0.0000000e+00\n".
    - Replaces ' 3"' with ' 2"'.
    - Replaces "XYZ" with "XY"."""
    x = open(mesh, "r").read()
    nb = x[: x.find("Topology")]
    rest = x[x.find("Topology") :]
    nb = nb.replace("0.0000000e+00\n", "\n")
    nb = nb.replace(' 3"', ' 2"')
    nb = nb.replace("XYZ", "XY")
    open(mesh, "w").write(nb + rest)


def vtp2xdmf(vtp):
    """
    Convert a VTP file to XDMF format and translate it to 2D.

    Parameters:
    vtp (str): The file path of the input VTP file. The file must have a .vtp extension.

    Returns:
    str: The file path of the converted XDMF file.

    Notes:
    - The function checks if the XDMF file already exists. If it does, it skips the conversion.
    - The function reads the VTP file, triangulates the mesh, sets the z-coordinates to 0 to make it 2D,
      and saves the mesh in XDMF format.
    - The function calls `conv3d_2d` to further process the XDMF file.
    - Prints messages indicating the status of the conversion.
    """
    assert vtp.endswith(".vtp")
    xd = vtp.replace(".vtp", ".xdmf")
    if os.path.exists(xd):
        logger.info(f"ANBA mesh {xd} already exists - skipping")
        return xd

    mesh = pv.read(vtp)
    tri = mesh.triangulate()  # fenics doesn't do mixed element types
    tri.save(vtp.replace(".vtp", "_tri.vtp"))
    tri.points[:, 2] = 0
    pv.save_meshio(xd, tri, data_format="XML")
    conv3d_2d(xd)
    logger.info(f"converted {vtp} to {xd}")
    return xd


def anba4_prep(section_meshes, parallel=True):
    """
    Convert a list of section meshes to XDMF format.

    This function takes a list of section meshes and converts each mesh to
    XDMF format using multiprocessing for parallel processing.

    Args:
        section_meshes (list): A list of section meshes to be converted.

    Returns:
        list: A list of converted section meshes in XDMF format.
    """
    logger.info("converting section meshes to XDMF")
    if parallel:
        p = multiprocessing.Pool()
        xds = p.map(vtp2xdmf, section_meshes)
        p.close()
        p.join()
    else:
        xds = []
        for i in section_meshes:
            xds.append(vtp2xdmf(i))
    logger.info("done converting section meshes to XDMF")
    return xds


def main():
    """
    Main function to translate section meshes from VTK to XDMF and 2D format.

    This function sets up an argument parser to handle command-line arguments
    for translating section meshes. It expects section meshes in the .vtp format
    and processes them accordingly.

    Command-line Arguments:
    sections (list of str): Section meshes in .vtp format.

    Returns:
    The result of the anba4_prep function, which processes the provided section meshes.
    """
    p = argparse.ArgumentParser(
        description=__doc__
        # "translate section meshes from vtk to XDMF and 2D"
    )
    p.add_argument("sections", nargs="*", help="section meshes in .vtp format")
    args = p.parse_args()
    return anba4_prep(args.sections)


if __name__ == "__main__":
    main()
