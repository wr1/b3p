# Blade modeling
## Building a blade model
A b3p blade model is built using the following command:
```bash 
b3p build blade_test.yml
```
This is subdivided into the following steps:
```bash
Usage: b3p COMMAND

╭─ Commands ───────────────────────────────────────────────────────────────────╮
│ 2d         Create 2D meshes for calculation of 6x6 matrices.                 │
│ build      Build the blade model: geometry, mesh, drape, and mass.           │
│ ccx        Run Calculix on this model                                        │
│ clean      Clean the working directory.                                      │
│ --help,-h  Display this message and exit.                                    │
│ --version  Display application version.                                      │
╰──────────────────────────────────────────────────────────────────────────────╯
```





# Input file format
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



