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

