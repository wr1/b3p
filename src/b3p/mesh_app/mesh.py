# Mesh building functions for mesh_app.

import os
import logging
import pickle
from copy import deepcopy as dc
import numpy as np
from pathlib import Path
import pyvista as pv
import vtk
from ..geometry.blade import blade
from ..geometry.blade_section import section
from ..geometry.loft_utils import load, interp, optspace
from ..geometry.splining import intp_c
from ..mesh.webs import write_web, build_webs
from ..geometry.geometry_blade_shape import blade_shape
from ..geometry.geometry_section import section as geometry_section
from ..geometry.geometry_web import web

logger = logging.getLogger(__name__)


def build_blade_mesh(config, workdir):
    """Build the blade mesh including webs."""
    pln = config["planform"]
    radii = np.linspace(0, 100, 100)
    web_inputs = config["mesh"]["webs"]
    base_vtp = workdir / "blade_geometry.vtp"
    web_intersections = build_webs(str(base_vtp), web_inputs, prefix="blade")
    prefix = "blade"
    pckfile = workdir / "blade_geometry.pck"
    outfile = workdir / "blade_mesh.vtp"
    build_mesh(
        str(pckfile),
        radii,
        web_inputs,
        web_intersections,
        prefix,
        outfile=str(outfile),
    )
    # Export web planes as JSON in workdir
    for web_name in web_intersections:
        data = {
            "name": web_name,
            "data": web_intersections[web_name],
            "points": {"lwp": [], "wwp": []},  # Points are not directly available here, but can be added if needed
        }
        with open(workdir / f"{web_name}.json", "w") as f:
            import json
            json.dump(data, f, indent=4)
    logger.info(f"Web planes exported to {workdir}")


def build_mesh(
    pckfile,
    radii,
    web_inputs,
    web_intersections,
    prefix,
    n_web_points=10,
    n_ch_points=120,
    outfile="out.vtp",
    added_datums=None,
    panel_mesh_scale=None,
):
    """Build the 3D mesh with webs linked to the shell."""
    if added_datums is None:
        added_datums = {}
    if panel_mesh_scale is None:
        panel_mesh_scale = []
    sections = pickle.load(open(pckfile, "rb"))
    weblist = [
        web(
            points=web_intersections[i],
            web_root=web_inputs[i]["z_start"],
            web_tip=web_inputs[i]["z_end"],
            web_name=f"{prefix}_{i}",
            coordinate=i,
            flip_normal=(web_inputs[i]["origin"][1] > 0),
        )
        for i in web_inputs
    ]
    nsec = []
    z = [i[0][2] for i in sections]
    for i in sections:
        r = i[0][2] - min(z)
        r_rel = (r - min(z)) / (max(z) - min(z))
        sec = geometry_section(r, r_rel, i, open_te=False)
        nsec.append(sec)
    blade = blade_shape(
        nsec,
        section_resolution=200,
        web_resolution=n_web_points,
        added_datums=added_datums,
        prefix=prefix,
    )
    for i in weblist:
        blade.set_web(i)
    blade.build_interpolated_sections(radii=radii, interpolation_type=2)
    blade.mesh(n_ch_points, panel_mesh_scale=panel_mesh_scale)
    blade.write_mesh(outfile)
    logger.info(f"Wrote blade mesh to {outfile}")
    return blade