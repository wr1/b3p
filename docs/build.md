# Blade modeling
## Building a blade model
A b3p blade model is built using the following command:
```bash 
b3p build blade_test.yml
```

This is subdivided into the following steps:
```bash
Usage: b3p build COMMAND [ARGS] [OPTIONS]

Build the blade model, geometry, mesh, drape
╭─ Commands ──────────────────────────────────────────────────╮
│ drape      Drape plies onto mesh                            │
│ geometry   Build blade geometry                             │
│ mesh       Mesh blade                                       │
│ --help,-h  Display this message and exit.                   │
│ --version  Display application version.                     │
╰─────────────────────────────────────────────────────────────╯
╭─ Parameters ────────────────────────────────────────────────╮
│ *  YML,--yml  b3p input file [required]                     │
╰─────────────────────────────────────────────────────────────╯
```
## Input file format
The top level keys in the input file are:
```yaml
general:    # run parameters
planform:   # blade planform
aero:       # airfoils and polars
mesh:       # meshing parameters including shear web position
mesh2d:     # 2D cross section meshing parameters
materials:  # material database
laminates:  # laminate plan
loads:      # loadcases
```

### general
```yaml
  prefix: test_blade    # prefix for output files
  workdir: temp_blade   # working directory
```


