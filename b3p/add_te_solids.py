#! /usr/bin/env python
import fire
import pyvista as pv
import pandas as pd


def load_vtu_file(file_path):
    # Load the VTU file
    mesh = pv.read(file_path)

    nr = 4
    df = pd.DataFrame(mesh.points, columns=["x", "y", "z"])
    df["d_te"] = mesh.point_data["d_te"]
    df["is_web"] = mesh.point_data["is_web"]
    df["d_rel_dist_from_te"] = mesh.point_data["d_rel_dist_from_te"]

    shell_pts = df[df.is_web == 0]

    # for g in shell_pts.groupby("z"):
    #     print(g[1].shape)

    grz = [g for g in shell_pts.groupby("z")]

    # grz = shell_pts.groupby("z")
    for cg, ng in zip(grz, grz[1:]):
        # Access current group
        print("Current Group:")
        print(cg[1].shape)

    # mesh.append_cells

    #     # print(current_group)

    #     # Access next group
    #     print("Next Group:")
    #     # print(next_group)

    # Perform operations on current and next group
    # ...

    # Example: Calculate the average d_te for current group
    # avg_d_te_current = current_group["d_te"].mean()
    # print("Average d_te for current group:", avg_d_te_current)

    # # Example: Calculate the average d_te for next group
    # avg_d_te_next = next_group["d_te"].mean()
    # print("Average d_te for next group:", avg_d_te_next)

    # print(shell_pts)
    # n = 100
    # for i in mesh.points[::n]:
    #     print(i)

    # Perform further operations on the mesh if needed
    # ...

    # Print some information about the loaded mesh
    print(mesh)


if __name__ == "__main__":
    fire.Fire(load_vtu_file)
