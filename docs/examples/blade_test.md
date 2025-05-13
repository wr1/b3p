# Running the `blade_test.yml` Example

This guide explains how to run the `blade_test.yml` example included in the **B3P** repository. The example demonstrates a complete workflow for processing a wind turbine blade model, including building the 3D geometry, performing 2D and 3D structural analyses, and generating output files.

## Prerequisites

Before running the example, ensure:

- B3P is installed (see [Installation](installation.md)).
- Dependencies (CalculiX, ANBA4, Conda environment `anba4-env`) are set up.
- The `examples` folder contains `blade_test.yml` and linked files (`materials.yml`, `loads.yml`, `laminates.yml`).
- Sufficient disk space (~500 MB) for output files.

## Running the Example

The `blade_test.yml` example processes a wind turbine blade model through geometry creation, meshing, 2D sectional analysis, and 3D finite element analysis (FEA). You can run it using the `testrun.py` script or individual B3P commands.

### Run all steps

```bash
cd examples
bash runex.sh
```


### Manual Execution

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

3. **Run 2D analysis**:
   ```bash
   b3p 2d blade_test.yml
   ```
   Generates 2D sectional meshes and runs ANBA4 analysis at specified radial positions.

4. **Run 3D FEA**:
   ```bash
   b3p ccx blade_test.yml
   ```
   Prepares and solves the 3D FEA using CalculiX for load cases (e.g., forward flap, edge).

5. **Run CCBlade analysis** (optional):
   ```bash
   b3p ccblade blade_test.yml
   ```
   Performs aerodynamic analysis.

## Generated Files

The workflow generates files in `temp_blade_portable`:

   - **Mesh Files** (`temp_blade_portable/mesh/`):
      - Airfoil data (e.g., `xs_test_blade_t_*.dat`).
      - 3D meshes (e.g., `test_blade.vtp`, `test_blade_shell.vtp`).
      - Web meshes (e.g., `test_blade_w0.vtp`).
      - Configuration files (e.g., `test_blade_variables.json`).

   - **Drape Files** (`temp_blade_portable/drape/`):
      - Material data (e.g., `material_map.json`, `ply_bom.csv`).
      - Draped meshes (e.g., `test_blade_shell_dr.vtu`).
      - Mass summaries (e.g., `test_blade_mass.csv`).
      - 2D sectional meshes (e.g., `msec_*.vtp`, `msec_*.xdmf`).

   - **FEA Files** (`temp_blade_portable/fea/`):
      - CalculiX inputs (e.g., `test_blade_ccx_mp_lc_forward_flap.inp`).
      - Results (e.g., `test_blade_ccx_mp_lc_forward_flap.vtu`).

   - **Aerodynamic Results**:
      - `ccblade_output.csv`, `ccblade_output.png`, `polars.png`.

## Workflow Summary

1. **Cleaning**: Clears previous output.
2. **Building**: Creates 3D geometry and meshes from YAML inputs.
3. **2D Analysis**: Generates and analyzes 2D sectional meshes.
4. **3D FEA**: Solves structural analysis for specified load cases.
5. **Aerodynamic Analysis**: Computes performance metrics with CCBlade.

## Notes
- The example generates ~200 files, including large meshes and result files.
- Typical outputs include a blade mass of ~57,321 kg and volume of ~66.38 m³.
- 2D analysis covers ~13 radial positions with convergence residuals (e.g., 1.77e-05).
- 3D FEA uses ~11,682 elements and solves in ~70 seconds per load case.
- Visualize results with Paraview or similar tools.

For troubleshooting, verify file paths, YAML syntax, and dependency installations. See [Usage](usage.md) for more details.

