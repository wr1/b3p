#! /usr/bin/env python

from b3p import yml_portable

import argparse
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import fmin
import numpy as np
from ccblade.ccblade import CCAirfoil, CCBlade
import math
import logging

logger = logging.getLogger(__name__)


def load_polar(pname):
    """
    Load a polar with a given name, interpolate to set alpha range.

    Args:
        pname (str): Path to the polar file.

    Returns:
        list: A list containing alpha, cl, cd, and cm values.
    """
    logger.info("loading polar", pname)
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


def plot_interpolated_polars(t, data, of="polars.png"):
    """
    Plot interpolated polars.

    Args:
        t (array-like): Thickness values.
        data (array-like): Polar data.
        of (str, optional): Output filename. Defaults to "polars.png".
    """
    fig, ax = plt.subplots(3, 1, figsize=(12, 19))
    for n, i in enumerate(t):
        alpha, cl, cd, cm = data[:, :, n]
        ax[0].plot(alpha, cl, label=f"t={t[n]:.3f}")
        ax[1].plot(alpha, cd)
        ax[2].plot(cd, cl)
        ax[0].set_ylabel("Cl")
        ax[1].set_ylabel("Cd")
        ax[2].set_ylabel("Cm")
        ax[0].legend()
    fig.savefig(of)


def plot_polars(polars, of="polars_in.png"):
    """
    Plot polars.

    Args:
        polars (list): List of polar data.
        of (str, optional): Output filename. Defaults to "polars_in.png".
    """
    fig, ax = plt.subplots(3, 1, figsize=(12, 16))
    for i in polars:
        alpha, cl, cd, cm = i[1]  # polars[i]
        ax[0].plot(alpha, cl, label=i[0])
        ax[1].plot(alpha, cd)
        ax[2].plot(cd, cl)
        ax[0].set_ylabel("Cl")
        ax[1].set_ylabel("Cd")
        ax[2].set_ylabel("Cl - Cd")
        ax[0].legend()
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
    fig.tight_layout()
    fig.savefig(of)


def interpolate_polars(polars, tnew, of=None):
    """
    Interpolate polars.

    Args:
        polars (list): List of polar data.
        tnew (pandas.Series): New thickness values.
        of (str, optional): Output filename. Defaults to None.

    Returns:
        list: List of interpolated CCAirfoil objects.
    """
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
        plot_polars(polars, of.replace(".png", "_in.png"))
        plot_interpolated_polars(tnew, data, of)

    return output_polars


class RotorOptimizer:
    """
    Optimize rotor performance.
    """

    def __init__(self, rotor, uinf, rated_power=20e7, omega=None, maxiter=5):
        """
        Initialize RotorOptimizer.

        Args:
            rotor (CCBlade): CCBlade rotor object.
            uinf (array-like): Inflow wind speeds.
            rated_power (float, optional): Rated power. Defaults to 20e7.
            omega (float, optional): Rotor speed. Defaults to None.
            maxiter (int, optional): Maximum iterations for optimization. Defaults to 5.
        """
        self.rotor = rotor
        self.uinf = uinf
        self.omega = omega
        self.rated_power = rated_power
        self.cache = {}
        self.latest = None
        self.maxiter = maxiter

    def evaluate(self, x, coefficients=False):
        """
        Evaluate rotor performance.

        Args:
            x (array-like): Optimization variables (omega, pitch).
            coefficients (bool, optional): Whether to return coefficients. Defaults to False.

        Returns:
            tuple: Power and evaluation outputs.
        """
        if self.omega == None:
            omega, pitch = x
        else:
            omega, pitch = self.omega, x[0]

        if (omega, pitch) not in self.cache:
            outputs, _ = self.rotor.evaluate(
                self.uinf, omega, pitch, coefficients=coefficients
            )
            self.cache[(omega, pitch)] = outputs

        P = np.mean(self.cache[(omega, pitch)]["P"])
        return P, self.cache[(omega, pitch)]

    def objective(self, x):
        """
        Objective function for optimization.

        Args:
            x (array-like): Optimization variables (omega, pitch).

        Returns:
            float: Absolute difference between rated power and evaluated power.
        """
        ev, _ = self.evaluate(x)
        return np.fabs(self.rated_power - ev)

    def optimize(self, initial_guess):
        """
        Optimize rotor performance.

        Args:
            initial_guess (array-like): Initial guess for optimization variables.

        Returns:
            array-like: Optimized values.
        """
        result = fmin(self.objective, initial_guess, maxiter=self.maxiter)

        rr, rdet = self.evaluate(result)
        return result


def tsr2omega(tsr, uinf, radius, max_tip_speed=95.0):
    """
    Convert tip speed ratio to rotor speed.

    Args:
        tsr (float): Tip speed ratio.
        uinf (float): Inflow wind speed.
        radius (float): Rotor radius.
        max_tip_speed (float, optional): Maximum tip speed. Defaults to 95.0.

    Returns:
        float: Rotor speed in RPM.
    """
    ts = np.minimum(uinf * tsr, max_tip_speed)
    return (ts * 60.0) / (2.0 * np.pi * radius)


def omega2tsr(omega, uinf, radius):
    """
    Convert rotor speed to tip speed ratio.

    Args:
        omega (float): Rotor speed in RPM.
        uinf (float): Inflow wind speed.
        radius (float): Rotor radius.

    Returns:
        float: Tip speed ratio.
    """
    return omega * 2.0 * np.pi * radius / (uinf * 60.0)


def plot_grid(num_plots, figsize=(15, 15)):
    """
    Create a grid of subplots.

    Args:
        num_plots (int): Number of subplots.
        figsize (tuple, optional): Figure size. Defaults to (15, 15).

    Returns:
        tuple: Figure and array of axes.
    """
    # Determine grid dimensions (try to get as square as possible)
    grid_size = math.isqrt(num_plots)

    # Adjust grid dimensions if necessary
    if grid_size * grid_size < num_plots:
        columns = grid_size
        rows = np.ceil(num_plots / grid_size).astype(int)
    else:
        rows = columns = grid_size

    fig, axs = plt.subplots(rows, columns, figsize=figsize)
    axs = axs.flatten()

    # Remove unused subplots
    for idx in range(len(axs) - 1, rows * columns):
        fig.delaxes(axs[idx])

    return fig, axs


def find_closest_x(x_values, evaluations, target, order):
    """
    Find the closest x value to a target value using polynomial interpolation.

    Args:
        x_values (array-like): X values.
        evaluations (array-like): Y values.
        target (float): Target value.
        order (int): Order of the polynomial.

    Returns:
        float: Closest x value.
    """
    assert len(x_values) == len(evaluations)
    # fit a polynomial of given order
    poly_coeffs = np.polyfit(x_values, evaluations, order)
    # create a polynomial function from the coefficients
    poly = np.poly1d(poly_coeffs)
    # create a dense range of x values for interpolation
    x_dense = np.linspace(x_values[0], x_values[-1], 10000)
    # find the interpolated x value where the absolute difference to the target is smallest
    x_closest = x_dense[np.argmin(np.abs(poly(x_dense) - target))]
    return x_closest


def plot_bladeloads(r, data_dict, of="bladeloads.png"):
    """
    Plot blade loads.

    Args:
        r (array-like): Radial positions.
        data_dict (dict): Dictionary of blade loads.
        of (str, optional): Output filename. Defaults to "bladeloads.png".
    """
    fig, axs = plot_grid(len(data_dict), figsize=(15, 15))
    # Iterate over the dictionary and plot each array
    for idx, (name, array) in enumerate(data_dict.items()):
        axs[idx].plot(r, array)
        axs[idx].set_title(name)

    fig.tight_layout()
    fig.savefig(of)
    logger.info(f"Saved {of}")


class controloptimize:
    """
    Optimize rotor control.
    """

    def __init__(self, rotor, max_tipspeed, rtip, rating, uinf, workdir):
        """
        Initialize controloptimize.

        Args:
            rotor (CCBlade): CCBlade rotor object.
            max_tipspeed (float): Maximum tip speed.
            rtip (float): Rotor tip radius.
            rating (float): Rated power.
            uinf (array-like): Inflow wind speeds.
            workdir (str): Working directory.
        """
        self.rotor = rotor
        self.max_tipspeed = max_tipspeed
        self.rating = rating
        self.uinf = uinf
        self.rtip = rtip
        self.workdir = workdir

    def control_opt_below_rated(
        self, starting_uinf=6, starting_tsr=10, starting_pitch=0
    ):
        """
        Optimize the control below rated power, optimize for both tsr and pitch, do this only for one wind speed (default 6m/s)

        Args:
            starting_uinf (int, optional): Wind speed to optimize for. Defaults to 6.
            starting_tsr (int, optional): Initial guess for tsr. Defaults to 10.
            starting_pitch (int, optional): Initial guess for pitch. Defaults to 0.
        """
        omega = tsr2omega(
            starting_tsr, uinf=6, radius=self.rtip, max_tip_speed=self.max_tipspeed
        )
        optimizer = RotorOptimizer(self.rotor, [starting_uinf])
        initial_guess = np.array([omega, starting_pitch])
        optimal_values = optimizer.optimize(initial_guess)
        logger.info(f"initial guess {initial_guess} optimal values {optimal_values}")
        init_val = optimizer.evaluate(initial_guess)[0]
        opt_val, optt = optimizer.evaluate(optimal_values, coefficients=True)

        logger.info(f" {init_val} {opt_val}, improvement {opt_val / init_val}")
        self.optimal_tsr = omega2tsr(optimal_values[0], starting_uinf, self.rtip)
        self.fine_pitch = optimal_values[1]

        loads, _ = self.rotor.distributedAeroLoads(
            self.uinf, optimal_values[0], optimal_values[1], 0
        )
        plot_bladeloads(
            self.rotor.r, loads, of=os.path.join(self.workdir, "ccblade_bladeloads.png")
        )
        logger.info(f"optimal tsr {self.optimal_tsr} {self.fine_pitch}")

    def control_opt_above_rated(self):
        """
        Optimize the control above rated power, keep tsr fixed and find the pitch that gives the rated power.
        """
        # first run through the powercurve with fine_pitch all the way
        self.omega = tsr2omega(
            self.optimal_tsr,
            uinf=self.uinf,
            radius=self.rtip,
            max_tip_speed=self.max_tipspeed,
        )
        self.tsr = omega2tsr(self.omega, uinf=self.uinf, radius=self.rtip)
        self.pitch = self.fine_pitch * np.ones_like(self.uinf)

        init_pc, _ = self.rotor.evaluate(
            self.uinf,
            self.omega,
            self.pitch,
            coefficients=False,
        )
        rotorplot(init_pc, self.uinf, of=os.path.join(self.workdir, "ccblade_init.png"))

        # find the wind speeds that are over rated power
        overrated = np.where(init_pc["P"] > self.rating)
        logger.info(f"overrated {overrated}, {self.uinf[overrated]}")
        upost = self.uinf[overrated]
        ompost = self.omega[overrated]
        pitch_over_rated = []
        closest_pitch = self.fine_pitch

        # keep tsr fixed and find the pitch that gives the rated power
        for i in np.array(list(zip(upost, ompost))):
            pitch = np.linspace(closest_pitch, closest_pitch + 15.0, 8)
            ui = np.ones_like(pitch) * i[0]
            oi = np.ones_like(pitch) * i[1]
            pg, _ = self.rotor.evaluate(ui, oi, pitch, coefficients=False)
            closest_pitch = find_closest_x(pitch, pg["P"], self.rating, 3)
            pitch_over_rated.append(closest_pitch)

        self.pitch[overrated] = pitch_over_rated

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

        # plot the result
        rotorplot(
            out_pc,
            self.uinf,
            labels=["P", "CP", "Mb", "T", "omega", "pitch", "tsr"],
            of=os.path.join(self.workdir, "ccblade_out.png"),
        )
        logger.info(f"pitch {self.pitch}")
        return out_pc


def rotorplot(op, uinf, labels=["P", "CP", "T", "Mb"], of="__temp.png"):
    """
    Plot rotor performance.

    Args:
        op (dict): Rotor performance data.
        uinf (array-like): Inflow wind speeds.
        labels (list, optional): List of labels to plot. Defaults to ["P", "CP", "T", "Mb"].
        of (str, optional): Output filename. Defaults to "__temp.png".
    """
    lab = [i for i in labels if i in op]

    fig, ax = plot_grid(len(lab), figsize=(10, 10))
    for n, i in enumerate(lab):
        ax[n].plot(uinf, op[i], label=f"{i} max={op[i].max():.2f}")
        ax[n].legend()
        ax[n].grid()
        ax[n].set_ylabel(i + "(W)" if i == "P" else "")
        ax[n].set_xlabel("uinf [m/s]")

    fig.tight_layout()
    fig.savefig(of)
    logger.info(f"saved {of}")


class ccblade_run:
    """
    Run CCBlade analysis.
    """

    def __init__(self, blade):
        """
        Initialize ccblade_run.

        Args:
            blade (str): Path to the blade YAML file.
        """
        self.dct = yml_portable.yaml_make_portable(blade)
        workdir = os.path.join(self.dct["general"]["workdir"], "mesh")
        bem = {
            "rated_power": 10e6,
            "polars": {},
            "B": 3,
            "rho": 1.225,
            "tilt": 5,
            "precone": 3,
            "shearExp": 0.1,
            "hubHt": 140.0,
            "mu": 1.81206e-5,
            "yaw": 0,
            "max_tipspeed": 95.0,
            "uinf": [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 20],
        }
        missing_keys = [key for key in bem if key not in self.dct["aero"]["bem"]]

        logger.info(
            f"the following keys are not specified in aero/bem: {missing_keys}, using defaults {bem}"
        )
        bem |= self.dct["aero"]["bem"]

        self.prefix = os.path.join(workdir, self.dct["general"]["prefix"])

        if "polars" not in bem:
            exit("no polars in blade file")

        if not os.path.isfile(f"{self.prefix}.pck"):
            exit("blade not built yet, run b3p build <blade.yml> first")

        plf = pd.read_csv(f"{self.prefix}_sca_50.csv", sep=";")

        plrs = sorted(
            [(i[0], load_polar(i[1])) for i in bem["polars"].items()],
            reverse=True,
        )
        iplr = interpolate_polars(
            plrs, plf.relative_thickness, of=os.path.join(workdir, "polars.png")
        )
        rhub, rtip = plf.z.iloc[0], plf.z.iloc[-1]
        self.rotor = CCBlade(
            plf.z - plf.z[0],
            plf.chord,
            plf.twist,
            iplr,
            rhub,
            rtip,
            B=bem["B"],
            rho=bem["rho"],
            mu=bem["mu"],
            precone=bem["precone"],
            tilt=bem["tilt"],
            yaw=bem["yaw"],
            shearExp=bem["shearExp"],
            hubHt=bem["hubHt"],
            derivatives=False,
        )

        logger.info(f"Rotor from {rhub} to {rtip}")

        self.copt = controloptimize(
            self.rotor,
            bem["max_tipspeed"],
            rtip,
            bem["rated_power"],
            uinf=np.array(bem["uinf"]),
            workdir=workdir,
        )

    def run(self):
        """
        Run the CCBlade analysis.
        """
        self.copt.control_opt_below_rated()
        output = self.copt.control_opt_above_rated()

        del output["W"]

        # Save the output to a CSV file
        workdir = self.dct["general"]["workdir"]
        df = pd.DataFrame(output).dropna()
        df.to_csv(os.path.join(workdir, "ccblade_output.csv"), sep=";")


def main():
    """
    Main function.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("blade", help="blade file")
    args = parser.parse_args()
    ccblade_run(args.blade).run()


if __name__ == "__main__":
    main()
