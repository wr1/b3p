"""Loadcase generation."""

import numpy as np
import logging

logger = logging.getLogger(__name__)


def get_loadcases(mesh, multiplier=1.0, buckling=False):
    """Get loadcases."""
    loadcases = {}

    for i in mesh.point_data:
        if i.startswith("lc_"):
            logger.info(f"loadcase {i}")
            multiplier = 1.0  # TODO fix for quadratic meshes
            if buckling:
                lbuf = f"** {i}\n*step\n*buckle\n5\n*cload\n"
            else:
                lbuf = f"** {i}\n*step\n*static\n*cload\n"

            ld = mesh.point_data[i] * multiplier
            for n, j in enumerate(ld):
                if j[0] ** 2 > 1e-8:
                    lbuf += "%i,1,%f\n" % (n + 1, j[0])
                if j[1] ** 2 > 1e-8:
                    lbuf += "%i,2,%f\n" % (n + 1, j[1])

            lbuf += "*node output,output=3d\nU\n*element output\nE,S\n*node print,nset=root,totals=yes\nrf\n*end step\n"

            loadcases[i] = lbuf

    return loadcases


def root_clamp(mesh):
    """Apply root clamp."""
    root = np.where(mesh.points[:, 2] == mesh.points[:, 2].min())
    lbuf = "*nset, nset=root\n"
    for n, j in enumerate(root[0]):
        lbuf += "%i" % (j + 1)
        if n % 16 == 0:
            lbuf += "\n"
        else:
            lbuf += ","
    lbuf += "\n"
    lbuf += "*boundary,op=new\nroot,1,3\n"
    return lbuf
