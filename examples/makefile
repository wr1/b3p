
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
	b3p_drape temp/plybook.pck temp/test_blade_web.vtp --key shell  --out temp/test_blade_shell1.vtu 
	b3p_drape temp/plybook.pck temp/test_blade_w0.vtp --key web0  --out temp/test_blade_w01.vtu 

combine: ## add the meshes together 
	b3p_meshcombine temp/test_blade_shell1.vtu  temp/test_blade_w01.vtu --out temp/joinedup_mesh.vtu

mesh2d: ## extract 2D meshes 
	b3p_2dmesh --vtu temp/joinedup_mesh.vtu --r "[2,10,20,30,40,50,60,70,80,90]" --var temp/test_blade.var

loads: ## assign nodal forces to the joined mesh for use in FEA
	b3p_add_load blade_test.yml temp/joinedup_mesh.vtu

ccx: ## export a grid to ccx format, solve and postprocess (requires ccx to be installed and be in PATH)
	b3p_export_ccx temp/joinedup_mesh.vtu --out temp/ccxmodel.inp	
	ccx temp/ccxmodel_lc_forward_flap
	ccx temp/ccxmodel_lc_forward_edge
	b3p_frd2vtu temp/ccxmodel_lc_forward_edge.frd --out temp/ccx_fe.vtu
	b3p_frd2vtu temp/ccxmodel_lc_forward_flap.frd --out temp/ccx_ff.vtu

anbasolve: ## solve using anba4
	python3 ../b3p/anba4_prep.py temp/msec_*.vtp
	python3 ../b3p/anba4_solve.py temp/msec_*.xdmf temp/material_map.json --debug

all: geom mesh plybook drape combine mesh2d loads ccx ## run all targets in sequence
