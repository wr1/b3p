# 2D cross section analysis

## Using ANBA to compute cross section properties (6x6 matrices)

```bash
# after running b3p build
# b3p build blade_test.yml
b3p 2d blade_test.yml
```

## Section analysis substeps

```bash
❯ b3p 2d -h
usage: b3p 2d [-h] [-r ROTZ] [-P] [-e ANBA_ENV] yml {mesh2d,run-anba4,clean} ...

positional arguments:
  yml                   Path to YAML config file
  {mesh2d,run-anba4,clean}
    mesh2d              Create 2D meshes
    run-anba4           Run ANBA4 on 2D meshes
    clean               Remove msec* files

options:
  -h, --help            show this help message and exit
  -r ROTZ, --rotz ROTZ  Rotation around Z-axis (degrees)
  -P, --no-parallel     Disable parallel processing for mesh2d
  -e ANBA_ENV, --anba-env ANBA_ENV
                        Conda environment for ANBA4
```

