#! /usr/bin/env python
import glob
import fire
import os
import pyvista as pv


class plot_ccx:
    def __init__(self, wdir):
        self.meshes = glob.glob(wdir + "/*ccx*.vtu")

    def plot3d(self):
        for i in self.meshes:
            self.__plot3d(pv.read(i), output_path=i.replace(".vtu", "_3d.png"))
        return self

    def __plot3d(self, mesh, output_path):

        # Deform the mesh using the ChannelDisplacement filter
        deformed_mesh = mesh.warp_by_vector("disp", factor=1)

        # Color the mesh using the third component of point data strain
        deformed_mesh.point_data["strain"] = mesh.point_data["strain"][:, 2]
        deformed_mesh.set_active_scalars("strain")

        deformed_mesh.rotate_y(-110, inplace=True)

        deformed_mesh.rotate_z(0, inplace=True)

        # Create a plot in isometric view
        plotter = pv.Plotter(off_screen=True, window_size=[1920, 1080])
        plotter.add_mesh(
            deformed_mesh,
            cmap="coolwarm",
            clim=[-5e-3, 5e-3],
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
        plotter.screenshot(
            output_path,
            transparent_background=True,
            return_img=False,
        )
        print(f"Saved {output_path}")

    def __str__(self):
        return ""

    def tab2d(self, num_bins=50):
        for i in self.meshes:
            self.__tab2d(
                pv.read(i), output_path=i.replace(".vtu", "_2d.txt"), num_bins=num_bins
            )
        return self

    def plot2d(self, num_bins=50):
        for i in self.meshes:
            self.__plot2d(
                pv.read(i), output_path=i.replace(".vtu", "_2d.png"), num_bins=num_bins
            )
        return self

    def __plot2d(self, mesh, output_path, num_bins=50):
        """"""
        z, ymin, ymax = digitize_strain_distribution(
            mesh.points[:, 2], mesh.point_data["strain"][:, 2], num_bins=num_bins
        )
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.plot(z, ymin, color="blue", label="compression")
        ax.plot(z, ymax, color="red", label="tension")

        ax.set_xlabel("z-coordinate")
        ax.set_ylabel("strain")
        ax.legend()
        ax.set_ylim(-7e-3, 7e-3)
        ax.grid(True)
        fig.savefig(output_path)
        print(f"Saved {output_path}")


if __name__ == "__main__":
    fire.Fire(plot_ccx)  # main()
