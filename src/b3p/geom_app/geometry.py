# Geometry building functions for geom_app.

from ..geometry.blade import blade
from ..geometry.loft_utils import optspace
from ..geometry.geometry_section import section as GeometrySection
from copy import deepcopy as dc
import os
import numpy as np
import pyvista as pv


def build_blade_geometry(config, workdir):
    """Perform blade 3d model generation based on b3p dictionary."""
    pln = config["planform"]
    blade_obj = blade(
        pln["chord"],
        pln["thickness"],
        pln["twist"],
        dc(pln["dx"]),
        pln["dy"],
        pln["z"],
        config["aero"]["airfoils"],
        chordwise_sampling=optspace(config["planform"]["npchord"]),
        np_spanwise=config["planform"]["npspan"],
    )

    # Use fixed file names in workdir
    blade_obj.mesh(str(workdir / "blade_geometry.vtp"))

    # Load the mesh and add full datums
    mesh = pv.read(workdir / "blade_geometry.vtp")
    n_points = blade_obj.np_chordwise
    for idx, sec in enumerate(blade_obj.sections):
        r = blade_obj.z[1][idx]  # Radius
        r_rel = idx / (len(blade_obj.sections) - 1)
        points = [[x, y, r] for x, y in zip(sec.x, sec.y)]
        gs = GeometrySection(r, r_rel, points, open_te=False)
        _, datums = gs.respline(n_points, [], {}, [])
        start_idx = idx * n_points
        for key, arr in datums.items():
            if key not in mesh.point_data:
                mesh.point_data[key] = np.zeros(len(mesh.points))
            mesh.point_data[key][start_idx : start_idx + n_points] = arr
    mesh.save(workdir / "blade_geometry.vtp")

    blade_obj.dump(str(workdir / "blade_geometry.pck"))
    blade_obj.export_variables(str(workdir / "blade_geometry_variables.json"))
    blade_obj.export_xfoil(
        prefix=os.path.join(workdir, "airfoil_out"),
    )
    blade_obj.plot(fname=str(workdir / "blade_geometry.png"))
    n_sections = 50
    blade_obj.to_table(
        np.linspace(0, 1, n_sections), str(workdir / "blade_geometry_sca_50.csv")
    )

    return blade_obj
