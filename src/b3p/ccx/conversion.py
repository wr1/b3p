"""FRD to VTU conversion utilities."""

import os
import glob
import logging
import frd2vtu

logger = logging.getLogger(__name__)


def has_later_vtu(frd):
    """Check if VTU is later than FRD."""
    output = frd.replace(".frd", ".vtu")
    if not os.path.exists(output):
        return False
    frd_time = os.path.getmtime(frd)
    vtu_time = os.path.getmtime(output)
    return frd_time < vtu_time


def all_frd2vtu(prefix):
    """Convert all FRD to VTU."""
    for i in glob.glob(f"{prefix}*.frd"):
        logger.info(f"processing {i}")
        if not has_later_vtu(i):
            frd2vtu(i)
        else:
            logger.info(f"skipping {i}")
