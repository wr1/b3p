# b3p (blade3 preprocessor)
Command line tools to create models for composite (wind turbine) blades. 
Note that this code does *not* include an FEA solver or any other solvers at this stage. 

## Functionality
- Interpolating a blade through a series of airfoils
- Aligning airfoils to create maximum building height over a straight shearweb
- Create rectangle grid for the blade
- Assign scalable *slab based* laminate plans 
- Assemble a shell mesh with laminate properties
- Write mesh information to VTK files

## How to run
```sh
pip install -e . 
cd examples
# build the geometry of the blade (closed geometry, triangles, no shearwebs)
make geom 
# create a structural mesh of the blade (including shearwebs)
make mesh 
# create a ply-based plybook out of slab-based input
make plybook 
# assign the plies to the mesh
make drape 
# combine the laminated shearwebs and shell meshes into a single grid with shared nodes
make combine 
```
## Output

![test_blade](https://user-images.githubusercontent.com/8971152/148471383-7f652a84-447a-4db0-81e2-2e27b1785745.png)

![3dblade_nplies](https://user-images.githubusercontent.com/8971152/148471469-61fb3efb-1789-4667-97b4-11b9e36d2e73.png)

## Input file format
TODO

