
# https://betterprogramming.pub/how-to-create-a-self-documenting-makefile-533ebf8f82e2
.PHONY: help
help: ## Show this help
	@egrep -h '\s##\s' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'


geom: ## create blade shell geometry (no shearwebs)
	b3p_blade_geometry blade_test.yml

mesh: ## create blade structure 
	b3p_blade_structure blade_test.yml

plybook: ## create a plybook out of slab definition
	b3p_plybook blade_test.yml --out temp/plybook.pck

drape: ## drape the shell and the shear web
	b3p_drape temp/plybook.pck temp/test_blade_web.vtp --key shell  --out temp/test_blade_shell.vtu 
	b3p_drape temp/plybook.pck temp/test_blade_w0.vtp --key web0  --out temp/test_blade_w0.vtu 

combine: ## add the meshes together 
	b3p_meshcombine temp/test_blade_shell.vtu  temp/test_blade_w0.vtu --out temp/joinedup_mesh.vtu

mesh2d: ## extract 2D meshes 
	b3p_2dmesh --vtu temp/joinedup_mesh.vtu --r "[5,20,35,50,65,80,90]" --var temp/test_blade.var

loads: ## assign nodal forces to the joined mesh for use in FEA
	b3p_add_load blade_test.yml temp/joinedup_mesh.vtu

ccx: ## export a grid to ccx format, solve and postprocess (requires ccx to be installed and be in PATH)
	b3p_export_ccx temp/joinedup_mesh.vtu --out temp/ccxmodel.inp	
	ccx temp/ccxmodel 
	b3p_frd2vtu temp/ccxmodel.frd 


all: geom mesh plybook drape combine mesh2d loads ccx 
