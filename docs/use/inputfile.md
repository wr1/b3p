# Input File Format

**B3P** uses YAML files to define wind turbine blade models, including geometry, materials, laminates, and loads. This guide describes the structure and key sections of a B3P input file.

## File Structure

A B3P input file is a YAML document with the following top-level sections:

- `general`: General settings (e.g., working directory, prefix).
- `aero`: Aerodynamic properties and airfoil data.
- `mesh`: Blade geometry and meshing parameters.
- `materials`: Material properties or path to material database.
- `laminates`: Composite laminate definitions.
- `loads`: Load cases for structural analysis.

## Example Input File

Below is a simplified example of a B3P input file:

```yaml
general:
  prefix: my_blade
  workdir: ./output_portable

aero:
  airfoils:
    0.18: path/to/airfoil_18.dat
    0.21: path/to/airfoil_21.dat
  bem:
    rated_power: 10e6  # 10 MW
    B: 3  # Number of blades
    rho: 1.225  # Air density
    max_tipspeed: 95.0
    uinf: [3, 5, 7, 9, 10, 11, 12, 13, 16, 20]

mesh:
  radii: [0.1, 5, 10, 20, 50, 100]
  n_web_points: 10
  n_chordwise_points: 120
  webs:
    web1:
      origin: [0, 0.1, 0]
      z_start: 0.1
      z_follow_blade: 50
      z_end: 100
      orientation: [0, 1, 0]
  bondline:
    material: adhesive
    width: [[0, 0], [0.5, 0.5], [1, 0.1]]

materials: path/to/materials.yml

laminates:
  slabs:
    shell:
      material: carbon_ud
      ply_thickness: 0.001
      slab: [[0, 0.01], [0.5, 0.02], [1, 0.015]]
      cover:
        chord: [0.1, 0.9, 0]
  datums:
    chord:
      xy: [[0, 0], [0.5, 0.5], [1, 1]]
      scalex: 1.0
      scaley: 1.0

loads:
  forward_flap:
    z: [0, 50, 100]
    mx: [0, 1e6, 2e6]
    my: [0, 0, 0]
    apply:
      radius: [0, 100]
```

## Section Details

### `general`
- `prefix`: Prefix for output files (e.g., `my_blade`).
- `workdir`: Working directory for output (e.g., `./output_portable`).

### `aero`
- `airfoils`: Dictionary mapping thickness ratios to airfoil data files (XFOIL format).
- `bem`: Blade element momentum (BEM) parameters for CCBlade:
  - `rated_power`: Turbine power rating (W).
  - `B`: Number of blades.
  - `rho`: Air density (kg/m³).
  - `max_tipspeed`: Maximum tip speed (m/s).
  - `uinf`: Array of wind speeds (m/s).

### `mesh`

- `radii`: List of radial positions for meshing (m).
- `n_web_points`: Number of points for web meshing.
- `n_chordwise_points`: Number of chordwise points for shell meshing.
- `webs`: Dictionary of web definitions:
  - `origin`: Web origin point [x, y, z].
  - `z_start`: Start of web along blade span (m).
  - `z_follow_blade`: Point where web follows blade geometry (m).
  - `z_end`: End of web (m).
  - `orientation`: Normal vector of web plane [x, y, z].
- `bondline`: Bondline configuration:
  - `material`: Material name for bondline.
  - `width`: List of [radius, width] pairs for bondline width variation.
- `coordinates` (optional): Additional datum points for meshing.

### `materials`

- Path to a YAML file defining material properties or a dictionary of materials.
- Example material file (`materials.yml`):
  ```yaml
  carbon_ud:
    name: Carbon UD
    e11: 135e9
    e22: 10e9
    e33: 10e9
    g12: 5e9
    g13: 5e9
    g23: 4e9
    nu12: 0.3
    nu13: 0.3
    nu23: 0.4
    rho: 1600
  adhesive:
    name: Adhesive
    E: 3e9
    nu: 0.35
    rho: 1200
  ```

### `laminates`

- `slabs`: Dictionary of laminate slabs:
  - `material`: Material name (must match `materials`).
  - `ply_thickness`: Thickness of each ply (m).
  - `slab`: List of [radius, thickness] pairs for thickness distribution.
  - `cover`: Chordwise coverage [start, end, increment].
  - `draping`: Either `plies` (individual plies) or `blocks` (core blocks).
- `datums`: Coordinate systems for coverage calculations.

### `loads`
- Dictionary of load cases (e.g., `forward_flap`):
  - `z`: Radial positions for load application (m).
  - `mx`, `my`: Moments about x and y axes (N·m).
  - `apply`: Conditions for load application (e.g., `radius: [min, max]`).

## Notes
- File paths (e.g., airfoils, materials) are relative to the YAML file unless absolute.
- Use `b3p yml_portable blade.yml` to create a portable version with embedded linked files.
- Validate YAML syntax using a linter or the B3P CLI.

For a complete example, see [Blade Test Example](examples/blade_test.md).

