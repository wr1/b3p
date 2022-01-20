# b3p (blade3 preprocessor)
Command line tools to create models for composite (wind turbine) blades. 

## Functionality
- Interpolating a blade through a series of airfoils
- Aligning airfoils to create maximum building height over a straight shearweb
- Creating a quad mesh for the blade
- Assigning scalable *slab based* laminate plans to the structure
- Assembling a shell mesh with laminate properties
- Writing mesh information to VTK files

## How to run
```sh
pip install -e . 
cd examples
make
```
which outputs the following targets
```
help                 Show this help
geom                 create blade shell geometry (no shearwebs)
mesh                 create blade structure 
plybook              create a plybook out of slab definition
drape                drape the shell and the shear web
combine              add the meshes together 
mesh2d               extract 2D meshes 
loads                assign nodal forces to the joined mesh for use in FEA
ccx                  export a grid to ccx format, solve and postprocess (requires ccx to be installed and be in PATH)
all                  run all targets in sequence
```
To run all steps type:
```
make all
``` 
## Output
Plot of geometric input parameters for the example blade
![test_blade](https://user-images.githubusercontent.com/8971152/148471383-7f652a84-447a-4db0-81e2-2e27b1785745.png)

Visualisation of the number of plies on the mesh using [Paraview](https://paraview.org)
![3dblade_nplies](https://user-images.githubusercontent.com/8971152/148471469-61fb3efb-1789-4667-97b4-11b9e36d2e73.png)

Visualisation of a 2d cross section mesh 
![mesh2d](https://user-images.githubusercontent.com/8971152/148645980-51c36e1a-89e1-469d-aeea-49bf5adf4070.png)

CalculiX results (very coarse mesh)
![anim3 0020](https://user-images.githubusercontent.com/8971152/150238270-3e057722-efca-4d07-bc28-27b786c456bd.png)


## Input file format
TODO


