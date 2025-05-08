# Running the `blade_test.yml` Test Run

This guide demonstrates how to execute the `blade_test.yml` example included in the **B3P** repository as a test run. This workflow processes a wind turbine blade model (SI122) through geometry creation, meshing, 2D and 3D structural analyses, and aerodynamic analysis, serving as a comprehensive test of B3P's functionality.

## Prerequisites

Before running the test, ensure:

- B3P is installed (see [Installation](../installation.md)).
- Dependencies (CalculiX, ANBA4, Conda environment `anba4-env`) are set up.
- The `examples` folder contains `blade_test.yml` and linked files (`materials.yml`, `loads.yml`, `laminates.yml`).
- Sufficient disk space (~500 MB) for output files.

## Running the Test

The `blade_test.yml` test run processes the SI122 blade model using a sequence of B3P commands. You can execute it using the `testrun.py` script or manually with CLI commands.

### Option 1: Using `testrun.py`

The `testrun.py` script automates the workflow:

```bash
cd examples
python testrun.py --example 1
```

This runs the full sequence of commands for the SI122 blade model.

### Option 2: Manual Execution

To run the workflow manually:

1. **Clean the workspace** (optional):
   ```bash
   b3p clean blade_test.yml
   ```
   Removes the `temp_blade_portable` directory.

2. **Build the 3D model**:
   ```bash
   b3p build blade_test.yml
   ```
   Constructs the blade geometry, loads airfoils (e.g., NACA 0017), and generates meshes.

3. **Drape composite plies**:
   ```bash
   b3p drape blade_test.yml
   ```
   Applies composite plies to the blade mesh.

4. **Run 2D analysis**:
   ```bash
   b3p 2d blade_test.yml
   ```
   Generates 2D sectional meshes and runs ANBA4 analysis at ~13 radial positions.

5. **Run 3D FEA**:
   ```bash
   b3p ccx blade_test.yml
   ```
   Prepares and solves 3D finite element analysis using CalculiX for load cases (e.g., forward flap, edge).

6. **Run CCBlade analysis**:
   ```bash
   b3p ccblade blade_test.yml
   ```
   Performs aerodynamic analysis.

## Expected Outputs

The test run generates files in the `temp_blade_portable` directory:

- **Mesh Files** (`temp_blade_portable/mesh/`):
  - Airfoil data (e.g., `xs_test_blade_t_*.dat`).
  - 3D meshes (e.g., `test_blade.vtp`, `test_blade_shell.vtp`).
  - Web meshes (e.g., `test_blade_w0.vtp`).
  - Configuration files (e.g., `test_blade_variables.json`).
- **Drape Files** (`temp_blade_portable/drape/`):
  - Material data (e.g., `material_map.json`, `ply_bom.csv`).
  - Draped meshes (e.g., `test_blade_shell_dr.vtu`).
  - Mass summaries (e.g., `test_blade_mass.csv`, ~57,321 kg).
  - 2D sectional meshes (e.g., `msec_*.vtp`, `msec_*.xdmf`).
- **FEA Files** (`temp_blade_portable/fea/`):
  - CalculiX inputs (e.g., `test_blade_ccx_mp_lc_forward_flap.inp`).
  - Results (e.g., `test_blade_ccx_mp_lc_forward_flap.vtu`, ~11,682 elements, ~70 seconds per load case).
- **Aerodynamic Results**:
  - `ccblade_output.csv`, `ccblade_output.png`, `polars.png`.

### Key Metrics

- **Blade Mass**: ~57,321 kg.
- **Blade Volume**: ~66.38 m³.
- **2D Analysis**: Covers ~13 radial positions with convergence residuals (e.g., 1.77e-05).
- **3D FEA**: Uses ~11,682 elements, solves in ~70 seconds per load case.

Visualize results using Paraview for `.vtp`, `.vtu`, and `.xdmf` files.

## Verifying the Test Run

To confirm the test run was successful:

- Check for the presence of output files in `temp_blade_portable` (~200 files).
- Verify key metrics (e.g., blade mass, convergence residuals) match expected values.
- Inspect log outputs for errors or warnings during execution.
- Visualize meshes and results in Paraview to ensure geometric and analytical accuracy.

## Troubleshooting

If the test run fails, consider:

- **File Not Found**: Verify `blade_test.yml` and linked files exist in the `examples` folder.
- **ANBA4 Errors**: Ensure the `anba4-env` Conda environment is activated:
  ```bash
  conda activate anba4-env
  ```
- **CalculiX Not Found**: Confirm `ccx` is in your PATH:
  ```bash
  which ccx
  ```
  Or specify the executable path with `--ccxexe`.
- **YAML Syntax Errors**: Validate `blade_test.yml` using a YAML linter or `b3p yml_portable blade_test.yml`.
- **Resource Issues**: Ensure sufficient disk space (~500 MB) and RAM for large meshes.

For persistent issues, check the [GitHub issues page](https://github.com/wr1/b3p/issues) or open a new issue with:

- Command output and error messages.
- Environment details (OS, Python version, B3P version, dependencies).
- Steps to reproduce the failure.

## Notes

- The test run generates ~200 files, including large meshes and result files.
- Use `b3p clean blade_test.yml` to remove temporary files after testing.
- The `testrun.py` script is recommended for consistency, as it ensures the correct sequence of commands.
- For detailed file formats and CLI options, see [Usage](../usage.md).

For a broader example context, refer to [Blade Test Example](blade_test.md).
