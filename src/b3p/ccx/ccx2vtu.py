#! /usr/bin/env python3

import pandas as pd
import pyvista
import numpy as np
import glob
import os
import re
import frd2vtu
import logging

logger = logging.getLogger(__name__)


def has_later_vtu(frd):
    output = frd.replace(".frd", ".vtu")
    if not os.path.exists(output):
        return False
    frd_time = os.path.getmtime(frd)
    vtu_time = os.path.getmtime(output)
    return frd_time < vtu_time


def all_frd2vtu(prefix):
    for i in glob.glob(f"{prefix}*.frd"):
        logger.info(f"processing {i}")
        if not has_later_vtu(i):
            frd2vtu(i)
        else:
            logger.info(f"skipping {i}")


def digitize_strain_distribution(z, strain, num_bins=100):
    # Check if strain has multiple components
    if strain.ndim == 1:
        strain = strain.reshape(-1, 1)
    elif strain.ndim != 2:
        raise ValueError("Invalid strain input. Expected 1D or 2D array.")

    num_components = strain.shape[1]

    # Bin the values in the first column into a specified number of bins
    zbin = np.linspace(z.min(), z.max(), num=num_bins + 1)
    bins = np.digitize(z, bins=zbin)

    # Find the maximum and minimum value in the second column for each bin and component
    max_vals = np.zeros((num_bins, num_components))
    min_vals = np.zeros((num_bins, num_components))

    for bin_num in range(1, num_bins + 1):
        for comp in range(num_components):
            max_vals[bin_num - 1, comp] = np.max(strain[bins == bin_num, comp])
            min_vals[bin_num - 1, comp] = np.min(strain[bins == bin_num, comp])

    return 0.5 * (zbin[:-1] + zbin[1:]), min_vals, max_vals


class ccx2vtu:
    def __init__(self, prefix, wildcard=""):
        self.prefix = prefix
        self.frds = glob.glob(f"{prefix}*lc*{wildcard}*.frd")

    def load_grids(self):
        logger.info(f"loading grids {self.frds}")
        self.grids = {}
        for i in self.frds:
            self.grids[i] = (
                pyvista.UnstructuredGrid(i.replace(".frd", ".vtu"))
                if has_later_vtu(i)
                else frd2vtu.frdbin2vtu(i)
            )

    def tabulate(self, n_bins=50):
        # Initialize an empty list to store the results
        if len(self.grids) == 0:
            logger.info(f"no grids found in {self.prefix}")
            return None

        # Loop through the file_list
        for grid in self.grids:
            # Read the input data from the current file
            input_data = self.grids[grid]
            # pd.read_csv(file_name)
            logger.info(grid)
            # logger.info(input_data.point_data.keys())

            strain0 = next(
                (
                    s
                    for s in input_data.point_data.keys()
                    if s.lower().find("strain") != -1
                ),
                None,
            )

            z = input_data.points[:, 2]
            strain = input_data.point_data[strain0]

            # Compute the digitized strain distribution for the current input
            zbin, min_vals, max_vals = digitize_strain_distribution(
                z, strain, num_bins=n_bins
            )

            # Combine min_vals and max_vals into a single array
            vals = np.hstack([min_vals, max_vals])

            # Create MultiColumns for the current input
            # logger.info(grid)
            columns = pd.MultiIndex.from_product(
                [
                    [re.search(r"_lc_(.+)\.frd", grid)[1]],
                    ["e_min", "e_max"],
                    ["xx", "yy", "zz", "xy", "yz", "xz"],
                ],
                names=["file", "type", "comp"],
            )

            # Create a DataFrame for the current input
            df = pd.DataFrame(vals, index=zbin, columns=columns)

            of = grid.replace(".frd", "_eps2d.pq")
            df.to_parquet(of)
            logger.info(f"written to {of}")

            # Append the current DataFrame to the results list
        #     results.append(df)

        # # Concatenate all the DataFrames in the results list
        # result_df = pd.concat(results)

        # return result_df
