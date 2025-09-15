# Optimizer classes for ccblade_app.

from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import fmin
from ccblade.ccblade import CCBlade
import os
import logging
from typing import List, Tuple, Dict, Optional, Any
from rich.live import Live
from rich.spinner import Spinner
from contextlib import redirect_stdout, redirect_stderr
import math

from .utils import tsr2omega, omega2tsr, find_closest_x
from .plots import (
    plot_interpolated_polars,
    plot_polars,
    plot_grid,
    plot_bladeloads,
    plot_moments,
    rotorplot,
)

logger = logging.getLogger(__name__)


class RotorOptimizer:
    """Optimize rotor performance."""

    def __init__(
        self,
        rotor: CCBlade,
        uinf: float,
        rated_power: float = 20e7,
        omega: Optional[float] = None,
        maxiter: int = 5,
    ):
        """Initialize the rotor optimizer with given parameters."""
        self.rotor = rotor
        self.uinf = uinf
        self.omega = omega
        self.rated_power = rated_power
        self.cache: Dict[tuple, Dict[str, Any]] = {}
        self.latest: Optional[Dict[str, Any]] = None
        self.maxiter = maxiter

    def evaluate(
        self, x: np.ndarray, coefficients: bool = False
    ) -> Tuple[float, Dict[str, Any]]:
        """Evaluate rotor performance for given variables."""
        if self.omega is None:
            omega, pitch = x
        else:
            omega, pitch = self.omega, x[0]

        if (omega, pitch) not in self.cache:
            with redirect_stdout(open(os.devnull, "w")), redirect_stderr(open(os.devnull, "w")):
                outputs, _ = self.rotor.evaluate(
                    self.uinf, omega, pitch, coefficients=coefficients
                )
            self.cache[(omega, pitch)] = outputs

        P = np.mean(self.cache[(omega, pitch)]["P"])
        return P, self.cache[(omega, pitch)]

    def objective(self, x: np.ndarray) -> float:
        """Compute objective function for optimization."""
        ev, _ = self.evaluate(x)
        return np.fabs(self.rated_power - ev)

    def optimize(self, initial_guess: np.ndarray) -> np.ndarray:
        """Optimize rotor variables using fmin."""
        result = fmin(self.objective, initial_guess, maxiter=self.maxiter)
        rr, rdet = self.evaluate(result)
        return result


class controloptimize:
    """Optimize rotor control settings."""

    def __init__(
        self,
        rotor: CCBlade,
        max_tipspeed: float,
        rtip: float,
        rating: float,
        uinf: np.ndarray,
        workdir: Path,
    ):
        """Initialize control optimizer with rotor parameters."""
        self.rotor = rotor
        self.max_tipspeed = max_tipspeed
        self.rating = rating
        self.uinf = uinf
        self.rtip = rtip
        self.workdir = workdir

    def control_opt_below_rated(
        self,
        starting_uinf: float = 6,
        starting_tsr: float = 10,
        starting_pitch: float = 0,
    ) -> None:
        """Optimize control for below rated power at a specific wind speed."""
        omega = tsr2omega(
            starting_tsr, uinf=6, radius=self.rtip, max_tip_speed=self.max_tipspeed
        )
        optimizer = RotorOptimizer(self.rotor, [starting_uinf])
        initial_guess = np.array([omega, starting_pitch])

        with Live(
            Spinner("dots", text="Optimize pitch and TSR below rated..."),
            refresh_per_second=10,
        ):
            with redirect_stdout(open(os.devnull, "w")), redirect_stderr(open(os.devnull, "w")):
                optimal_values = optimizer.optimize(initial_guess)
                init_val = optimizer.evaluate(initial_guess)[0]
                opt_val, optt = optimizer.evaluate(optimal_values, coefficients=True)

        logger.info(f" {init_val} {opt_val}, improvement {opt_val / init_val}")
        self.optimal_tsr = omega2tsr(optimal_values[0], starting_uinf, self.rtip)
        self.fine_pitch = optimal_values[1]
        logger.info(f"optimal tsr {self.optimal_tsr} {self.fine_pitch}")

    def control_opt_above_rated(self) -> Dict[str, Any]:
        """Optimize control for above rated power by adjusting pitch."""
        self.omega = tsr2omega(
            self.optimal_tsr,
            uinf=self.uinf,
            radius=self.rtip,
            max_tip_speed=self.max_tipspeed,
        )
        self.tsr = omega2tsr(self.omega, uinf=self.uinf, radius=self.rtip)
        self.pitch = self.fine_pitch * np.ones_like(self.uinf)

        with Live(
            Spinner("dots", text="Above rated initial guess..."),
            refresh_per_second=10,
        ):
            with redirect_stdout(open(os.devnull, "w")), redirect_stderr(open(os.devnull, "w")):
                init_pc, _ = self.rotor.evaluate(
                    self.uinf,
                    self.omega,
                    self.pitch,
                    coefficients=False,
                )
        rotorplot(init_pc, self.uinf, of=self.workdir.parent / "ccblade_init.png")
        overrated = np.where(init_pc["P"] > self.rating)
        logger.info(f"overrated {overrated}, {self.uinf[overrated]}")
        upost = self.uinf[overrated]
        ompost = self.omega[overrated]
        pitch_over_rated = []
        closest_pitch = self.fine_pitch
        with Live(
            Spinner(
                "dots", text="Evaluating for adjusted pitch to find fixed power..."
            ),
            refresh_per_second=10,
        ):
            for i in np.array(list(zip(upost, ompost))):
                pitch = np.linspace(closest_pitch, closest_pitch + 15.0, 8)
                ui = np.ones_like(pitch) * i[0]
                oi = np.ones_like(pitch) * i[1]
                with redirect_stdout(open(os.devnull, "w")), redirect_stderr(open(os.devnull, "w")):
                    pg, _ = self.rotor.evaluate(ui, oi, pitch, coefficients=False)
                closest_pitch = find_closest_x(pitch, pg["P"], self.rating, 3)
                pitch_over_rated.append(closest_pitch)
        self.pitch[overrated] = pitch_over_rated
        with Live(
            Spinner("dots", text="Final rotor evaluation..."), refresh_per_second=10
        ):
            with redirect_stdout(open(os.devnull, "w")), redirect_stderr(open(os.devnull, "w")):
                out_pc, _ = self.rotor.evaluate(
                    self.uinf,
                    self.omega,
                    self.pitch,
                    coefficients=True,
                )
        out_pc["omega"] = self.omega
        out_pc["tsr"] = self.tsr
        out_pc["pitch"] = self.pitch
        out_pc["uinf"] = self.uinf
        rotorplot(
            out_pc,
            self.uinf,
            labels=["P", "CP", "Mb", "T", "omega", "pitch", "tsr"],
            of=self.workdir.parent / "ccblade_out.png",
        )
        logger.info(f"pitch {self.pitch}")
        return out_pc

    def compute_bladeloads(self) -> None:
        """Compute and plot blade loads for all operating points."""
        loads_list = []
        uinf_list = []
        flapwise_moments = []
        edgewise_moments = []
        r = self.rotor.r
        for ui, om, pi in zip(self.uinf, self.omega, self.pitch):
            with redirect_stdout(open(os.devnull, "w")), redirect_stderr(open(os.devnull, "w")):
                loads, _ = self.rotor.distributedAeroLoads(ui, om, pi, 0)
            loads_list.append(loads)
            uinf_list.append(ui)
            # Compute moments
            Np = loads['Np']
            Tp = loads['Tp']
            moment_flap = np.trapz(Np * r, r)
            moment_edge = np.trapz(Tp * r, r)
            flapwise_moments.append(moment_flap)
            edgewise_moments.append(moment_edge)
        plot_bladeloads(self.rotor.r, loads_list, uinf_list, of=self.workdir.parent / "ccblade_bladeloads.png")
        # Save table output
        all_loads = {}
        for ui, loads in zip(uinf_list, loads_list):
            ui_str = f"{ui:.1f}"
            for key, arr in loads.items():
                if key in ['Np', 'Tp']:
                    col_name = f"{key}_{ui_str}"
                    all_loads[col_name] = arr
        df_all = pd.DataFrame(all_loads, index=self.rotor.r)
        df_all.index.name = 'r'
        df_all.to_csv(self.workdir.parent / "ccblade_bladeloads.csv")
        logger.info(f"Saved blade loads table")
        # Save moments
        combined_rms = np.sqrt(np.array(flapwise_moments)**2 + np.array(edgewise_moments)**2)
        moments_dict = {
            'uinf': uinf_list,
            'flapwise_moment': flapwise_moments,
            'edgewise_moment': edgewise_moments,
            'combined_rms': combined_rms
        }
        df_moments = pd.DataFrame(moments_dict)
        df_moments.to_csv(self.workdir.parent / "ccblade_moments.csv", index=False)
        logger.info(f"Saved moments table")
        # Plot moments
        plot_moments(self.rotor.r, loads_list, uinf_list, {'flapwise': np.array(flapwise_moments), 'edgewise': np.array(edgewise_moments), 'combined_rms': combined_rms}, of=self.workdir.parent / "ccblade_moments.png")
