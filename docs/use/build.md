# Blade modeling
## Building a blade model
A b3p blade model is built using the following command.
```bash 
b3p build blade_test.yml
```
## Substeps
The build command runs a series of substeps that form a 3D blade mesh with assigned materials. 
```bash
❯ b3p build -h
usage: b3p build [-h] [-n] yml {geometry,mesh,drape,mass,apply-loads} ...

positional arguments:
  yml                   Path to YAML config file
  {geometry,mesh,drape,mass,apply-loads}
    geometry            Build blade geometry
    mesh                Mesh blade structure
    drape               Drape plies onto mesh
    mass                Calculate blade mass
    apply-loads         Apply loads to mesh

options:
  -h, --help            show this help message and exit
  -n, --no-bondline     Exclude bondline
```

## Interpolation of planform parameters
The blade planform is defined in the `planform` subsection of the input file.
```yaml
planform:
  npchord: 100
  npspan: 70
  # offset in prebend direction, negative is upwind
  dx: [[0.0, 0.0], [0.2, 0.0], [0.4, 0.0], [1.0, -5.0]]
  # offset in chord direction
  dy: [[0.0, 0.0], [0.01, 0], [0.3, 0.0, 0.5860957549908372], [0.43, 0.5292475355828817],
    [0.56, 0.4981144528302349], [0.69, 0.43552499282773205],
    [0.82, 0.32879129475154056], [0.95, 0.3334681541732368]]
  z: [[0, 3], [1, 126]]
  chord: [[0.0, 5], [0.03, 5.06], [0.2, 6.5], [0.4, 4.4], [0.6, 2.8], [0.75, 2.0], [
      0.94, 1.4], [0.985, 0.8], [1, 0.1]]
  thickness: [[0.0, 1.0], [0.05, 0.97], [0.2, 0.53], [0.31, 0.39], [1.0, 0.17]]
  twist: [[0.0, 10], [0.2, 2], [0.4, 1], [0.6, 0], [0.8, -0.5], [0.9, -1], [1.0, 0]]
```
Which gets interpolated. 
![Test blade](../assets/images/test_blade.png)

## Creation of 3D models for the shell and webs



<!-- ## Creation of a 2D blade model -->




