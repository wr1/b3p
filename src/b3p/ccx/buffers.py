"""Node and element buffer generation."""

import numpy as np
import vtk
import logging

logger = logging.getLogger(__name__)


def nodebuffer(grid):
    """Generate node buffer."""
    nodes = np.column_stack((np.arange(1, len(grid.points) + 1), grid.points))
    logger.info(f"exporting {len(nodes)} nodes")
    return (
        "*node,nset=nall\n"
        + "\n".join([f"{int(n)},{x:f},{y:f},{z:f}" for n, x, y, z in nodes])
        + "\n"
    )


def element_buffer(grid):
    """Generate element buffer."""
    conn = grid.cells_dict
    # logger.info(f"{conn}")
    extypes = [23]

    vtk_ccx = {23: "s8r"}

    buf = ""
    for tp in extypes:
        ccxtype = vtk_ccx[tp]
        for n, i in enumerate(conn[tp]):
            conn[tp][n] = np.array(i) + 1
            buf += f"*element,type={ccxtype},elset=e{n + 1}\n"
            buf += f"{n + 1},{','.join(map(str, conn[tp][n]))}\n"

    return buf


def orientation_buffer(grid, add_centers=False):
    """Generate orientation buffer."""
    shells = grid.extract_cells(grid.cells_dict.get(23, []))
    buf = ""
    x_dirs = shells.cell_data["x_dir"]
    y_dirs = shells.cell_data["y_dir"]
    centers = shells.cell_data["centers"]
    num_cells = shells.GetNumberOfCells()

    for n in range(num_cells):
        xdir, ydir = x_dirs[n], y_dirs[n]
        center = centers[n]
        buf += "*orientation,name=or%i,system=rectangular\n" % (n + 1)
        if add_centers:
            coords = np.concatenate([xdir + center, ydir + center, center], axis=0)
            buf += ",".join(format(k, ".4g") for k in coords.tolist()) + "\n3,0\n"
        else:
            coords = np.concatenate([xdir, ydir], axis=0)
            buf += ",".join(format(k, ".4g") for k in coords.tolist()) + "\n"

    return buf
