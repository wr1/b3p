#! /usr/bin/env python3
import pyvista as pv
import sys
import numpy as np
import fire


def visualize_slab_thickness(vtu_file, output_image=None):
    if output_image is None:
        output_image = vtu_file.replace(".vtu", ".png")
    # Read the unstructured grid from the VTU file
    grid = pv.read(vtu_file)

    egrid = grid.extract_geometry()

    # Extract the slab arrays
    slab_arrays = [name for name in grid.cell_arrays.keys() if name.startswith("slab_")]

    if not slab_arrays:
        print("No arrays with names starting with 'slab_' were found in the VTU file.")
        sys.exit(1)

    # Create a plotter object
    plotter = pv.Plotter(off_screen=True)

    # Define the x-offset for the horizontal stacking
    x_offset = 1.3 * abs(grid.bounds[2] - grid.bounds[3])
    # 1.5 times the z-size of the original grid

    # Loop over the slab arrays and add them to the plotter with the appropriate offset and rotation
    for i, slab_name in enumerate(slab_arrays):
        slab_grid = egrid.copy(deep=True)
        slab_grid.cell_arrays[slab_name] = grid.cell_arrays[slab_name]
        # Set the active scalars to the current slab array
        slab_grid.rotate_y(90)  # Rotate the grid around the y-axis by 90 degrees

        if "web" in slab_name:
            non_zero_mask = slab_grid.cell_arrays[slab_name] != 0
            slab_grid = slab_grid.extract_cells(non_zero_mask)
            slab_grid.rotate_x(90)  # Rotate the grid around the y-axis by 90 degrees

        slab_grid.translate([0, i * x_offset, 0])  # Apply horizontal stacking

        plotter.add_mesh(
            slab_grid,
            show_edges=False,
            cmap="viridis",
            clim=None,
            scalars=slab_name,
            scalar_bar_args={
                "title": slab_name.split("thickness_")[-1],
                "position_x": 0.85,
                "position_y": 0.05 + (i * 0.9 / len(slab_arrays)),
                "width": 0.1,
                "height": 0.1,
                "color": "black",
                "n_labels": 2,
                "title_font_size": 22,
            },
        )

    # Set the camera position
    plotter.view_xy()
    plotter.camera.zoom(1.6)
    # Save the high-resolution image
    plotter.screenshot(
        output_image, transparent_background=True, window_size=[1920, 1080]
    )

    print(f"Saved the high-resolution image to: {output_image}")


if __name__ == "__main__":
    fire.Fire(visualize_slab_thickness)