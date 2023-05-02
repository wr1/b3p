#! /usr/bin/env python

from b3p import geometry_section
from b3p import geometry_blade_shape
from b3p import geometry_web
import pickle


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

    if added_datums is None:
        added_datums = {}
    if panel_mesh_scale is None:
        panel_mesh_scale = []
    sections = pickle.load(open(pckfile, "rb"))

    weblist = [
        geometry_web.web(
            points=web_intersections[i],
            web_root=web_inputs[i]["z_start"],
            web_tip=web_inputs[i]["z_end"],
            web_name=f"{prefix}_{i}.txt",
            coordinate=i,
            # normal=web_inputs[i]["normal"] if "normal" in web_inputs[i] else (0, 1, 0),
            flip_normal=(web_inputs[i]["origin"][1] > 0),
        )
        for i in web_inputs
    ]
    nsec = []
    z = [i[0][2] for i in sections]

    # loop over the sections and add them to a section list
    for i in sections:
        r = i[0][2] - min(z)
        r_rel = (r - min(z)) / (max(z) - min(z))
        sec = geometry_section.section(r, r_rel, i, open_te=False)
        nsec.append(sec)

    # create a blade loft from the sections
    blade = geometry_blade_shape.blade_shape(
        nsec,
        section_resolution=200,
        web_resolution=n_web_points,
        added_datums=added_datums,
    )

    for i in weblist:
        blade.set_web(i)

    blade.build_interpolated_sections(radii=radii, interpolation_type=2)

    # mesh with a given number of points around the circumference
    blade.mesh(n_ch_points, panel_mesh_scale=panel_mesh_scale)
    blade.write_mesh(outfile)
    print(f"** wrote mesh to {outfile}")

    return blade
