"""Material database handling for CCX."""

import json
import numpy as np
import logging

logger = logging.getLogger(__name__)


def material_db_to_ccx(materials, matmap=None, force_iso=False):
    """Find the material db and write properties to a ccx block"""
    mat_db = None
    if os.path.isfile(matmap):  # check if the material map file is there
        mm1 = json.load(open(matmap, "r"))
        mm = mm1["map"]
        mat_db = mm1["matdb"]
    else:
        exit("no material map defined")

    mm_inv = {v: k for k, v in mm.items()}

    matblock = ""
    for i in materials:
        if i > 1e-6:
            material_properties = mat_db[mm_inv[int(i)]]
            matblock += (
                f"** material: {mm_inv[int(i)]} {i} {material_properties['name']}\n"
            )

            if "C" in material_properties and not force_iso:
                logger.info(
                    f"{material_properties['name']} is assumed to be orthotropic"
                )
                C = np.array(material_properties["C"])
                matblock += "** orthotropic material\n"
                matblock += "*material,name=m%i\n*elastic,type=ortho\n" % i
                D = C
                D[0, 3] = C[0, 5]
                D[0, 5] = C[0, 3]
                D[1, 3] = C[1, 5]
                D[1, 5] = C[1, 3]
                D[2, 3] = C[2, 5]
                D[2, 5] = C[2, 3]
                D[3, 3] = C[5, 5]
                D[5, 5] = C[3, 3]
                matblock += (
                    f"{D[0, 0]:.4g},{D[0, 1]:.4g},{D[1, 1]:.4g},"
                    + f"{D[0, 2]:.4g},{D[1, 2]:.4g},{D[2, 2]:.4g},"
                    + f"{D[3, 3]:.4g},{D[4, 4]:.4g},\n"
                    + f"{D[5, 5]:.4g},293\n"
                )
            elif "Ex" in material_properties and not force_iso:
                logger.info(f"{material_properties['name']} has engineering constants")
                matblock += "** orthotropic material\n"
                matblock += (
                    "*material,name=m%i\n*elastic,type=engineering constants\n" % i
                )
                matblock += (
                    f"{material_properties['Ex']:.4g},{material_properties['Ey']:.4g},{material_properties['Ez']:.4g},"
                    + f"{material_properties['nuxy']:.4g},{material_properties['nuxz']:.4g},{material_properties['nuyz']:.4g},"
                    + f"{material_properties['Gxy']:.4g},{material_properties['Gxz']:.4g},\n"
                    + f"{material_properties['Gyz']:.4g},293\n"
                )
            else:
                logger.info(f"{material_properties['name']} is assumed to be isotropic")
                nu = min(
                    0.45,
                    max(
                        0.1,
                        (
                            float(material_properties["nu"])
                            if "nu" in material_properties
                            else material_properties["nuxy"]
                        ),
                    ),
                )
                E = float(
                    material_properties["Ex"]
                    if "Ex" in material_properties
                    else material_properties["E"]
                )
                matblock += "** isotropic material\n"
                matblock += "*material,name=m%i\n*elastic,type=iso\n" % i
                matblock += f"{E:.4g},{nu:.4g},293\n"

    return matblock
