"""Element set formatting and computation."""

import numpy as np
import logging

logger = logging.getLogger(__name__)


def format_eset(name, eids):
    """Format element set."""
    out = f"*elset,elset={name}\n"
    for i in range(len(eids)):
        out += f"{eids[i]}"
        out += "\n" if (i % 16 == 15) else ","
    if out[-1] == ",":
        out = out[:-1] + "\n"
    return out


def compute_ply_groups(grid, prefix):
    """Compute ply groups."""
    gr = ""
    n = 1
    for i in grid.cell_data:
        if i.startswith(prefix):
            eids = np.where(grid.cell_data[i][:, 1] > 0)[0] + 1
            gr += format_eset(i, eids)
            n += 1
    return gr


def compute_slab_groups(grid, prefix):
    """Compute slab groups."""
    gr = ""
    for i in grid.cell_data:
        if i.startswith(prefix):
            eids = np.where(grid.cell_data[i] > 0)[0] + 1
            gr += format_eset(i, eids)
    return gr
