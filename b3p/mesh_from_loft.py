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
    added_datums={},
    panel_mesh_scale=[],
):

    sections = pickle.load(open(pckfile, "rb"))

    # create webs, using parametric coordinates, the point lists are
    # (rel_coordinate_along_web,first_rel_chordwise_point,second_rel_chordwise_point)
    # web_root is the lowest r location of the web, web_tip is the highest

    weblist = []
    c = 0
    for i in web_inputs:
        weblist.append(
            geometry_web.web(
                points=web_intersections[i],
                web_root=web_inputs[i]["z_start"],
                web_tip=web_inputs[i]["z_end"],
                web_name="%s_%s.txt" % (prefix, i),
                coordinate=i,
                flip_normal=(web_inputs[i]["origin"][1] > 0),
            )
        )
        c += 1

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
    print("# writing to %s" % outfile)
    blade.write_mesh(outfile)
    print("# writing mesh done")
