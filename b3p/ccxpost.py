#! /usr/bin/env python
import glob
import fire
import os
import pyvista as pv
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from cycler import cycler


class plot_ccx:
    def __init__(self, wdir, wildcard=""):
        self.wdir = wdir
        self.meshes = glob.glob(wdir + f"/*ccx*{wildcard}*.vtu")

    def plot3d(self):
        for i in self.meshes:
            mesh = pv.read(i)
            self.__plot3d(mesh, output_path=i.replace(".vtu", "_3d.png"))
            del mesh
        return self

    def __plot3d(self, mesh, output_path):
        ts = np.unique([i.split("_")[1] for i in mesh.point_data.keys()])

        for i in ts:
            isdisp = i.startswith("1.00000")
            # Deform the mesh using the ChannelDisplacement filter
            deformed_mesh = mesh.warp_by_vector(
                f"disp_{ts[0]}", factor=1  # if isdisp else None
            )

            # Color the mesh using the third component of point data strain
            deformed_mesh.point_data["strain"] = mesh.point_data[f"strain_{i}"][:, 2]
            deformed_mesh.set_active_scalars("strain")

            deformed_mesh.rotate_y(-110, inplace=True)

            deformed_mesh.rotate_z(0, inplace=True)

            # Create a plot in isometric view
            plotter = pv.Plotter(off_screen=True, window_size=[1920, 1080])
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
            plotter.camera.zoom(2.0)

            # Save the plot to a high-resolution image
            of = output_path.replace(".png", f"_{i}.png")
            plotter.screenshot(
                of,
                transparent_background=True,
                return_img=False,
            )
            print(f"Saved {of}")
            # clear up memory
            plotter.close()
            del plotter

    def __str__(self):
        return ""

    def plot2d(self, wildcard=""):
        pqs = glob.glob(os.path.join(self.wdir, f"*{wildcard}*eps2d.pq"))

        df = pd.concat([pd.read_parquet(i) for i in pqs])
        fig, ax = plt.subplots(figsize=(12, 8))
        # colors = ["red", "blue", "green", "orange"]
        ax.set_prop_cycle(
            cycler(color=["red", "blue", "green", "orange"])
            + cycler(linestyle=["-", "--", ":", "-."])
        )
        for i in df:
            if i[2] == "zz":
                ax.plot(
                    df[i].index,
                    df[i].values,
                    label=" ".join(i[:2]),
                )
        ax.set_xlabel("z-coordinate")
        ax.set_ylabel("strain")
        ax.legend()
        ax.set_ylim(-7e-3, 7e-3)
        ax.grid(True)
        output_path = pqs[0].replace(".pq", ".png")
        fig.savefig(output_path, dpi=300)  # , transparent=True)
        print(f"Saved {output_path}")


if __name__ == "__main__":
    fire.Fire(plot_ccx)  # main()
