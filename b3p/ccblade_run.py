#! /usr/bin/env python

from utils import yml_portable
import fire
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import fmin
import numpy as np
from ccblade.ccblade import CCAirfoil, CCBlade
import math


def load_polar(pname):
    "load a polar with a given name, interpolate to set alpha range"
    print("loading polar", pname)
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
    fig, ax = plt.subplots(3, 1, figsize=(6, 8))
    for i in polars:
        alpha, cl, cd, cm = i[1]  # polars[i]
        ax[0].plot(alpha, cl, label=i[0])
        ax[1].plot(alpha, cd)
        ax[2].plot(cd, cl)
        ax[0].set_ylabel("Cl")
        ax[1].set_ylabel("Cd")
        ax[2].set_ylabel("Cl - Cd")
        ax[0].legend()
    fig.savefig(of)


def interpolate_polars(polars, tnew, of=None):
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
    def __init__(self, rotor, uinf, rated_power=20e7, omega=None, maxiter=5):
        self.rotor = rotor
        self.uinf = uinf
        self.omega = omega
        self.rated_power = rated_power
        self.cache = {}
        self.latest = None
        self.maxiter = maxiter

    def evaluate(self, x, coefficients=False):
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
        ev, _ = self.evaluate(x)
        return np.fabs(self.rated_power - ev)  # + np.sum(penalties)

    def optimize(self, initial_guess):
        result = fmin(self.objective, initial_guess, maxiter=self.maxiter)

        rr, rdet = self.evaluate(result)
        # print(f"{self.uinf} {self.omega} {self.rated_power} result {result}")
        # print(f"result {rr} {rdet}")
        return result


def tsr2omega(tsr, uinf, radius, max_tip_speed=95.0):
    ts = np.minimum(uinf * tsr, max_tip_speed)
    return (ts * 60.0) / (2.0 * np.pi * radius)


def omega2tsr(omega, uinf, radius):
    return omega * 2.0 * np.pi * radius / (uinf * 60.0)


def plot_grid(num_plots, figsize=(15, 15)):
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
        # print("deleting", idx)
        fig.delaxes(axs[idx])

    return fig, axs  # , rows, columns


def find_closest_x(x_values, evaluations, target, order):
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
    fig, axs = plot_grid(len(data_dict), figsize=(15, 15))
    # Iterate over the dictionary and plot each array
    for idx, (name, array) in enumerate(data_dict.items()):
        axs[idx].plot(r, array)
        axs[idx].set_title(name)

    fig.tight_layout()
    fig.savefig(of)
    print(f"Saved {of}")


class controloptimize:
    def __init__(self, rotor, max_tipspeed, rtip, rating, uinf, workdir):
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

        starting_uinf: wind speed to optimize for
        starting_tsr: initial guess for tsr
        starting_pitch: initial guess for pitch
        """
        omega = tsr2omega(
            starting_tsr, uinf=6, radius=self.rtip, max_tip_speed=self.max_tipspeed
        )
        optimizer = RotorOptimizer(self.rotor, [starting_uinf])
        initial_guess = np.array([omega, starting_pitch])
        optimal_values = optimizer.optimize(initial_guess)
        print(f"initial guess {initial_guess} optimal values {optimal_values}")
        init_val = optimizer.evaluate(initial_guess)[0]
        opt_val, optt = optimizer.evaluate(optimal_values, coefficients=True)

        print(f" {init_val} {opt_val}, improvement {opt_val/init_val}")
        self.optimal_tsr = omega2tsr(optimal_values[0], starting_uinf, self.rtip)
        self.fine_pitch = optimal_values[1]

        loads, _ = self.rotor.distributedAeroLoads(
            self.uinf, optimal_values[0], optimal_values[1], 0
        )
        plot_bladeloads(
            self.rotor.r, loads, of=os.path.join(self.workdir, "bladeloads.png")
        )
        print(f"optimal tsr {self.optimal_tsr} {self.fine_pitch}")

    def control_opt_above_rated(self):
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
        rotorplot(init_pc, self.uinf, of=os.path.join(self.workdir, "init.png"))

        # find the wind speeds that are over rated power
        overrated = np.where(init_pc["P"] > self.rating)
        print(f"overrated {overrated}, {self.uinf[overrated]}")
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

        # plot the result
        rotorplot(
            out_pc,
            self.uinf,
            labels=["P", "CP", "Mb", "T", "omega", "pitch", "tsr"],
            of=os.path.join(self.workdir, "out.png"),
        )
        print(out_pc.keys())
        print(f"pitch {self.pitch}")
        return out_pc


def rotorplot(op, uinf, labels=["P", "CP", "T", "Mb"], of="__temp.png"):
    lab = [i for i in labels if i in op]

    # fig, ax = plt.subplots(len(lab), 1, figsize=(10, 10))
    fig, ax = plot_grid(len(lab), figsize=(10, 10))
    # print(len(ax), len(lab), rows, columns)
    for n, i in enumerate(lab):
        ax[n].plot(uinf, op[i], label=f"{i} max={op[i].max():.2f}")
        ax[n].legend()
        ax[n].grid()
        ax[n].set_ylabel(i)

    fig.tight_layout()
    fig.savefig(of)
    print(f"saved {of}")


class ccblade_run:
    def __init__(self, blade):
        self.dct = yml_portable.yaml_make_portable(blade)
        workdir = self.dct["general"]["workdir"]
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
            "uinf": [3, 5, 7, 9, 10, 11, 12, 13, 16, 20],
        }
        missing_keys = [key for key in bem if key not in self.dct["aero"]["bem"]]

        print(
            f"the following keys are not specified in aero/bem: {missing_keys}, using defaults {bem}"
        )
        bem |= self.dct["aero"]["bem"]

        self.prefix = os.path.join(
            self.dct["general"]["workdir"], self.dct["general"]["prefix"]
        )

        if "polars" not in bem:
            exit("no polars in blade file")

        if not os.path.isfile(f"{self.prefix}.pck"):
            exit("blade not run yet")

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

        print(f"Rotor from {rhub} to {rtip}")

        copt = controloptimize(
            self.rotor,
            bem["max_tipspeed"],
            rtip,
            bem["rated_power"],
            uinf=np.array(bem["uinf"]),
            workdir=workdir,
        )

        copt.control_opt_below_rated()

        copt.control_opt_above_rated()

    def __str__(self):
        return ""


# def main(blade):
#     ccblade_run(blade=blade)


# def ccblade_run_test():
#     ccblade_run()


# if __name__ == "__main__":
#     fire.Fire(main)


def main():
    fire.core.Display = lambda lines, out: print(*lines, file=out)
    fire.Fire(ccblade_run)


if __name__ == "__main__":
    main()
