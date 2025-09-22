"""Convert mesh to CCX input."""

import pyvista as pv
import numpy as np
import vtk
import time
import logging
from .material_db import material_db_to_ccx
from .element_sets import compute_ply_groups, compute_slab_groups
from .buffers import nodebuffer, element_buffer, orientation_buffer
from .loadcases import get_loadcases, root_clamp
from .shell_sections import make_shell_section

logger = logging.getLogger(__name__)


def zero_midside_loads(mesh):
    """Zero midside loads."""
    if mesh.celltypes[0] == 23:
        conn = mesh.cell_connectivity.reshape(
            (
                mesh.GetNumberOfCells(),
                int(mesh.cell_connectivity.shape[0] / mesh.GetNumberOfCells()),
            )
        )
        midsides = conn[:, 4:8].flatten()
        for i in mesh.point_data:
            if i.startswith("lc_"):
                mesh.point_data[i][midsides] *= 0.0
    return mesh


def mesh2ccx(
    vtu,
    out="test.inp",
    matmap="temp/material_map.json",
    merge_adjacent_layers=True,
    zeroangle=False,
    single_step=False,
    quadratic=True,
    add_centers=False,
    force_isotropic=False,
    # export_hyperworks=False,
    export_plygroups=False,
    buckling=False,
    meshonly=False,
    bondline=False,  # Added to accept bondline argument
):
    """Convert VTU to CCX input file."""
    logger.info(f"Converting {vtu} to ccx input file {out}")
    grid = pv.read(vtu)
    gr = grid.threshold(value=(1e-6, 1e9), scalars="thickness")
    gr.cell_data["centers"] = gr.cell_centers().points

    logger.info(f"Exporting {gr.GetNumberOfCells()} elements")
    if quadratic:
        lf = vtk.vtkLinearToQuadraticCellsFilter()
        lf.SetInputData(gr)
        lf.Update()
        quad = lf.GetOutput()
        mesh = zero_midside_loads(pv.UnstructuredGrid(quad))
        mesh.save(vtu.replace(".vtu", "_quad.vtu"))
    else:
        mesh = pv.UnstructuredGrid(gr)

    buf = nodebuffer(mesh)

    if quadratic:
        buf += element_buffer(mesh)
    else:
        buf += element_buffer(mesh)

    if meshonly:
        of = out.replace(".inp", "_meshonly.inp")
        open(of, "w").write(buf)
        logger.info(f"written to {of}")
        return [of]

    buf += "*elset,elset=Eall,GENERATE\n%i,%i\n" % (1, mesh.GetNumberOfCells())

    if export_plygroups:
        plygroups = compute_ply_groups(mesh, "ply_")
        slabgroups = compute_slab_groups(mesh, "slab_thickness_")
        buf += plygroups + slabgroups

    plykeys = [i for i in mesh.cell_data if i.startswith("ply_")]
    plydat = np.stack([mesh.cell_data[i] for i in plykeys])
    materials = np.unique(plydat[:, :, 0])

    buf += orientation_buffer(mesh, add_centers)

    logger.info("made orientation buffer")
    matblock = material_db_to_ccx(materials, matmap=matmap, force_iso=force_isotropic)

    buf += "** START MATERIALS\n"
    buf += matblock
    buf += "** END MATERIALS\n"

    tic = time.perf_counter()

    blx = [
        make_shell_section(i, plydat[:, i, :], merge_adjacent_layers, zeroangle)
        for i in range(plydat.shape[1])
    ]
    logger.info("made shell sections")

    nplies = np.array([i[0] for i in blx])
    nplmax = nplies.max()
    npxid = np.where(nplies == nplmax)[0]
    logger.info(f"max number of plies: {nplmax}")
    # logger.info(f"associated stack \n{blx[npxid[0]][1]}")

    toc = time.perf_counter()
    logger.info(f"time spent creating shell sections {toc - tic:f}")
    logger.info(f"mesh {mesh}")
    comps = "".join(
        f"*shell section,composite,elset=e{n + 1},offset=-.5"
        + (f",orientation=or{n + 1}\n" if zeroangle else "\n")
        + i[1]
        for n, i in enumerate(blx)
    )
    buf += "** START SHELL SECTIONS\n"
    buf += comps
    buf += "** END SHELL SECTIONS\n"
    buf += root_clamp(mesh)

    loadcases = get_loadcases(mesh, buckling=buckling)

    logger.info(f"written loadcases {loadcases.keys()}")
    output_files = []
    if single_step:
        output = buf + "".join(loadcases.values())
        open(out, "w").write(output)
        logger.info(f"written ccx input file with all loadcases to {out}")
        output_files.append(out)
    else:
        for i in loadcases:
            output = buf + loadcases[i]
            of = out.replace(".inp", f"_{i}.inp")
            open(of, "w").write(output)
            output_files.append(of)
            logger.info(f"written ccx input file to {of}")

    return output_files
