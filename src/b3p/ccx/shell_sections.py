"""Shell section creation."""

import numpy as np
import logging

logger = logging.getLogger(__name__)


def make_shell_section(elem_id, plyarray, merge_adjacent_plies=True, zero_angle=True):
    """Make shell section."""
    plies = []
    filtered_plyarray = plyarray[plyarray[:, 1] > 1e-6]

    for j in filtered_plyarray:
        if plies and plies[-1][1] == j[0] and merge_adjacent_plies:
            plies[-1][0] += j[1] * 1e-3
        else:
            plies.append([j[1] * 1e-3, j[0]])

    if zero_angle:
        section_string = "".join("%f,,m%i,0\n" % tuple(i) for i in plies)
    else:
        section_string = "".join(
            "%f,,m%i,or%i\n" % tuple(i + [elem_id + 1]) for i in plies
        )

    return len(plies), section_string
