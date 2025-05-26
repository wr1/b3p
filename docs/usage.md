# Usage Guide

This guide explains how to use **B3P** (Blade Preprocessor) via its command-line interface (CLI) and Python API. B3P processes wind turbine blade models defined in YAML files, enabling geometry creation, meshing, and analyses.

## Command-Line Interface (CLI)

B3P provides a CLI with commands for building models, running analyses, and cleaning workspaces. The general syntax is:

```bash
b3p <command> <yaml_file> [options]
```

### Key Commands

1. **Build**: Constructs the 3D blade model and meshes.
   ```bash
   b3p build blade.yml
   ```
   Subcommands:
   - `geometry`: Build blade geometry only.
   - `mesh`: Generate meshes.
   - `drape`: Drape composite plies.
   - `mass`: Calculate blade mass.
   - `apply-loads`: Apply loads to the mesh.
   Example:
   ```bash
   b3p build blade.yml --no-bondline
   ```

2. **2D Analysis**: Generates 2D sectional meshes and runs ANBA4 analysis.
   ```bash
   b3p 2d blade.yml
   ```
   Subcommands:
   - `mesh2d`: Create 2D meshes.
   - `run-anba4`: Run ANBA4 analysis.
   - `clean`: Remove temporary `msec*` files.
   Options:
   - `--rotz <degrees>`: Rotate meshes around Z-axis.
   - `--anba-env <env>`: Specify Conda environment for ANBA4.
   Example:
   ```bash
   b3p 2d blade.yml --rotz 10 --anba-env anba4-env
   ```

3. **CCX (CalculiX)**: Runs 3D finite element analysis.
   ```bash
   b3p ccx blade.yml
   ```
   Subcommands:
   - `prep`: Prepare CalculiX input files.
   - `solve`: Solve the FEA problem.
   - `post`: Post-process results.
   - `plot`: Generate 2D/3D plots.
   Options:
   - `--bondline`: Include bondline meshes.
   - `--nproc <n>`: Number of processes for solving.
   - `--wildcard <pattern>`: Filter input/result files.
   Example:
   ```bash
   b3p ccx blade.yml --nproc 4 --bondline
   ```

4. **CCBlade**: Performs aerodynamic analysis.
   ```bash
   b3p ccblade blade.yml
   ```

5. **Clean**: Removes temporary files from the working directory.
   ```bash
   b3p clean blade.yml
   ```

### Example Workflow

To process a blade model from start to finish:

```bash
b3p clean blade.yml
b3p build blade.yml
b3p 2d blade.yml
b3p ccx blade.yml
b3p ccblade blade.yml
```

## Python API

B3P can be used programmatically for custom workflows. Below is an example of running a CCBlade analysis:

```python
from b3p.ccblade_run import ccblade_run

# Run CCBlade analysis
analysis = ccblade_run("blade.yml")

# Access results
print(analysis.results)
```

For advanced control optimization:

```python
from b3p.ccblade_run import ccblade_run, controloptimize

blade = ccblade_run("blade.yml")
rotor = blade._setup_rotor()
copt = controloptimize(
    rotor,
    max_tipspeed=95.0,
    rtip=blade.rtip,
    rating=15e6,
    uinf=np.array([3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]),
    workdir="./custom_output"
)
copt.control_opt_below_rated(starting_tsr=8, starting_pitch=2)
output = copt.control_opt_above_rated()
```

## Output Files

B3P generates files in the working directory (e.g., `output_portable`):
- **Mesh Files**: `.vtp`, `.vtu` for 3D geometry; `.xdmf` for 2D meshes.
- **Drape Files**: `material_map.json`, `ply_bom.csv`, `mass.csv` for composite data.
- **FEA Results**: `.frd`, `.vtu` for CalculiX outputs.
- **CCBlade Results**: `ccblade_output.csv`, `ccblade_output.png`, `polars.png`.

Use visualization tools like Paraview to inspect `.vtp`, `.vtu`, and `.xdmf` files.

## Troubleshooting

- **File not found**: Verify YAML file paths and linked files (e.g., airfoils, materials).
- **ANBA4 errors**: Ensure the `anba4-env` Conda environment is activated.
- **CalculiX failures**: Check `ccx` executable path and input file syntax.

Refer to [Examples](examples/blade_test.md) for detailed workflows.

