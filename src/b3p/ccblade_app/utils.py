# Utility functions for ccblade_app.

from pathlib import Path
import numpy as np
import pandas as pd
from ccblade.ccblade import CCAirfoil
import os
import logging
from typing import List, Tuple, Dict, Optional, Any

logger = logging.getLogger(__name__)


def load_polar(pname: Path) -> Tuple[List[float], np.ndarray, np.ndarray, np.ndarray]:
    """Load and interpolate a polar by name to a set alpha range."""
    logger.info(f"loading polar {pname}")
    if not os.path.isfile(pname):
        raise IOError(f"Polar {pname} not found")

    with open(pname, "r") as f:
        lines = f.readlines()

    read_start_index = next(
        i for i, line in enumerate(lines) if line.startswith("-180")
    )
    read_end_index = next((i for i, line in enumerate(lines) if "EOT" in line), None)

    lst = [
        [float(j) for j in line.split()]
        for line in lines[read_start_index:read_end_index]
    ]
    g = list(zip(*lst))
    alpha_new = (
        list(np.linspace(-180, -20, 25))
        + list(np.linspace(-19, 20, 29))
        + list(np.linspace(21, 180, 25))
    )
    cl = np.interp(alpha_new, g[0], g[1])
    cd = np.interp(alpha_new, g[0], g[2])
    cm = np.interp(alpha_new, g[0], g[3])
    return [alpha_new, cl, cd, cm]


def interpolate_polars(
    polars: List[tuple], tnew: np.ndarray, of: Optional[Path] = None
) -> List[CCAirfoil]:
    """Interpolate polars to new thickness values and return CCAirfoil objects."""
    indices = range(len(polars[0][1][0]))
    t = [polar[0] for polar in polars]
    data = np.array(
        [
            [
                np.interp(
                    np.flip(tnew.values),
                    np.flip(t),
                    np.flip([polar[1][i][idx] for polar in polars]),
                )
                for idx in indices
            ]
            for i in range(4)
        ]
    )
    output_polars = [
        CCAirfoil(
            Re=[1e6],
            alpha=data[0, :, i],
            cl=data[1, :, i],
            cd=data[2, :, i],
            cm=data[3, :, i],
        )
        for i in reversed(range(len(tnew)))
    ]
    if of:
        from .plots import plot_polars, plot_interpolated_polars
        plot_polars(polars, of=of.with_name(of.stem + "_in" + of.suffix))
        plot_interpolated_polars(tnew, data, of=of)
    return output_polars


def tsr2omega(
    tsr: float, uinf: float, radius: float, max_tip_speed: float = 95.0
) -> float:
    """Convert tip speed ratio to rotor speed in RPM."""
    ts = np.minimum(uinf * tsr, max_tip_speed)
    return (ts * 60.0) / (2.0 * np.pi * radius)


def omega2tsr(omega: float, uinf: float, radius: float) -> float:
    """Convert rotor speed to tip speed ratio."""
    return omega * 2.0 * np.pi * radius / (uinf * 60.0)


def find_closest_x(
    x_values: np.ndarray, evaluations: np.ndarray, target: float, order: int
) -> float:
    """Find closest x value to target using polynomial interpolation."""
    assert len(x_values) == len(evaluations)
    poly_coeffs = np.polyfit(x_values, evaluations, order)
    poly = np.poly1d(poly_coeffs)
    x_dense = np.linspace(x_values[0], x_values[-1], 10000)
    x_closest = x_dense[np.argmin(np.abs(poly(x_dense) - target))]
    return x_closest
