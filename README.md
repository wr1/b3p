[![Deploy](https://github.com/wr1/b3p/actions/workflows/publish.yml/badge.svg)](https://github.com/wr1/b3p/actions/workflows/publish.yml)[![Test](https://github.com/wr1/b3p/actions/workflows/test.yml/badge.svg)](https://github.com/wr1/b3p/actions/workflows/test.yml)![PyPI](https://img.shields.io/pypi/v/b3p)

# b3p 
Command line tools to create models for composite (wind turbine) blades. 

## Functionality
- Building 3D models of wind turbine blades
- Creating a quad mesh for the blade
- Assigning scalable *slab based* laminate plans to the structure
- Assembling a shell mesh with laminate properties
- Writing mesh information to VTK files

## Install
From pypi
```sh
pip install b3p
``` 
From source
```sh
git clone https://github.com/wr1/b3p.git
cd b3p 
pip install -e .
``` 
## How to run
```sh 
usage: b3p [-h] {build,ccx,2d,ccblade,clean} ...

Blade Design CLI

positional arguments:
  {build,ccx,2d,ccblade,clean}
    build               Build the full blade model
    ccx                 Run Calculix operations
    2d                  2D mesh and ANBA4 operations
    ccblade             Run CCBlade analysis
    clean               Clean working directory

options:
  -h, --help            show this help message and exit
```
In a cloned repository:
```sh
cd examples
# build the blade model, geometry, mesh, drape
b3p build blade_test.yml
# run the ccx fea analysis
b3p ccx blade_test.yml  
```

## Output
Plot of geometric input parameters for the example blade
![test_blade](https://user-images.githubusercontent.com/8971152/148471383-7f652a84-447a-4db0-81e2-2e27b1785745.png)

Visualisation of the number of plies on the mesh using [Paraview](https://paraview.org)
![3dblade_nplies](https://user-images.githubusercontent.com/8971152/148471469-61fb3efb-1789-4667-97b4-11b9e36d2e73.png)

Visualisation of a 2d cross section mesh 
![mesh2d](https://user-images.githubusercontent.com/8971152/148645980-51c36e1a-89e1-469d-aeea-49bf5adf4070.png)

[CalculiX](http://www.dhondt.de/) results (very coarse mesh)
![zstrain](https://user-images.githubusercontent.com/8971152/151350188-0a6f31bf-5f0e-457b-b6cb-438bb10b4c91.png)


