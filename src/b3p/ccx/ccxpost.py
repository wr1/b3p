#! /usr/bin/env python
import glob
import pyvista as pv
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from cycler import cycler
import os
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def get_quadrant(vector):
    ang = np.degrees(np.arctan2(vector[1], vector[0]))

    if np.isclose(ang, -90, atol=15):
        return -180
    elif np.isclose(ang, 90, atol=15):
        return 0
    elif np.isclose(ang, 0, atol=15):
        return 90
    elif np.isclose(ang, 180, atol=15):
        return -90
    else:
        return 45


class plot_ccx:
    def __init__(self, wdir: Path, wildcard: str = ""):
        self.wdir = wdir
        self.meshes = [
            i for i in wdir.glob((f"fea/*ccx_lc*{wildcard}*.vtu".replace("**", "*")))
        ]
        logger.info(self.meshes)

    def plot3d(self):
        for i in self.meshes:
            mesh = pv.read(i)

            self.__plot3d(mesh, output_path=str(i).replace(".vtu", "_3d.png"))
            del mesh
        return self

    def __plot3d(self, mesh, output_path):
        ts = [
            i.split("_")[1]
            for i in mesh.point_data.keys()
            if i.lower().startswith("disp")
        ]
        DISPKEY = "DISP_" + ts[0] if ts else None
        STRAINKEY = "STRAIN_" + ts[0] if ts else None

        logger.info(f"Found time steps: {ts}")

        for i in ts:
            isdisp = i.startswith("0.00000")
            # Deform the mesh using the ChannelDisplacement filter
            deform = mesh.point_data[f"DISP_{ts[0]}"]
            deformed_mesh = mesh.warp_by_vector(
                f"DISP_{ts[0]}",
                factor=1,  # if isdisp else None
            )

            if isdisp:  # f"strain_{i}" in mesh.point_data.keys():
                # Color the mesh using the third component of point data strain
                deformed_mesh.point_data["strain"] = mesh.point_data[f"STRAIN_{i}"][
                    :, 2
                ]
                deformed_mesh.set_active_scalars("strain")
            else:
                # Color the mesh using the third component of point data strain
                deformed_mesh.point_data["disp"] = mesh.point_data[f"DISP_{i}"].sum(
                    axis=1
                )
                deformed_mesh.set_active_scalars("disp")

            nrm = np.linalg.norm(deform, axis=1)

            whr = np.where(nrm == nrm.max())
            quadrant = get_quadrant(deform[whr[0][0]][:2])

            deformed_mesh.rotate_z(quadrant, inplace=True)

            deformed_mesh.rotate_y(-110, inplace=True)

            # Create a plot in isometric view
            plotter = pv.Plotter(off_screen=True, window_size=[1920, 600])
            plotter.add_mesh(
                deformed_mesh,
                cmap="coolwarm",
                clim=[-5e-3, 5e-3] if isdisp else None,
                scalar_bar_args={
                    "position_x": 0.1,
                    "position_y": 0.05,
                    "width": 0.8,
                    "height": 0.06,
                    "color": "black",
                    "n_labels": 2,
                    "title_font_size": 22,
                },
            )
            plotter.view_isometric()
            plotter.camera.parallel_projection = True
            plotter.camera.zoom(3.5)

            # Save the plot to a high-resolution image
            of = output_path.replace(".png", f"_{i}.png")
            plotter.screenshot(
                of,
                transparent_background=True,
                return_img=False,
            )
            logger.info(f"Saved {of}")
            # clear up memory
            plotter.close()
            del plotter

    def __str__(self):
        return ""

    def plot2d(self, wildcard=""):
        pqs = glob.glob(os.path.join(self.wdir, "fea", f"*{wildcard}*eps2d.pq"))

        logger.info(f"Found {len(pqs)} parquet files for 2D plotting: {pqs}")

        for pq in pqs:
            df = pd.read_parquet(pq)
            fig, ax = plt.subplots(figsize=(12, 8))
            ax.set_prop_cycle(
                cycler(color=["red", "blue", "green", "orange"])
                + cycler(linestyle=["-", "--", ":", "-."])
            )
            for i in df:
                if i[2] == "zz":
                    ax.plot(
                        df[i].index.values,
                        df[i].values,
                        label=" ".join(i[:2]),
                    )
            ax.set_xlabel("z-coordinate")
            ax.set_ylabel("strain")
            ax.legend()
            # ax.set_ylim(-7e-3, 7e-3)
            ax.grid(True)
            output_path = pq.replace(".pq", ".png")
            fig.savefig(output_path, dpi=300)  # , transparent=True)
            logger.info(f"Saved {output_path}")
            plt.close(fig)
