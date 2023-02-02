
# https://betterprogramming.pub/how-to-create-a-self-documenting-makefile-533ebf8f82e2
.PHONY: help
help: ## Show this help
	@egrep -h '\s##\s' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'


YML = blade_test_portable.yml
TMP= temp_portable

geom: ## create blade shell geometry (no shearwebs)
	b3p_blade_geometry $(YML) 

mesh: ## create blade structure 
	b3p_blade_structure $(YML)

plybook: ## create a plybook out of slab definition
	b3p_plybook $(YML) --out $(TMP)/plybook.pck

drape: ## drape the shell and the shear web
	b3p_drape $(TMP)/plybook.pck $(TMP)/test_blade_web.vtp --key shell  --out $(TMP)/test_blade_shell1.vtu 
	b3p_drape $(TMP)/plybook.pck $(TMP)/test_blade_w4.vtp --key web4  --out $(TMP)/test_blade_w04.vtu 
	b3p_drape $(TMP)/plybook.pck $(TMP)/test_blade_w3.vtp --key web3  --out $(TMP)/test_blade_w03.vtu 

combine: ## add the meshes together 
	b3p_meshcombine $(TMP)/test_blade_shell1.vtu   $(TMP)/test_blade_w03.vtu --out $(TMP)/joinedup_mesh.vtu

# $(TMP)/test_blade_w04.vtu
mesh2d: ## extract 2D meshes 
	b3p_2dmesh --vtu $(TMP)/joinedup_mesh.vtu --r "[2,10,20,30,40,50,60,70,80,90]" --var $(TMP)/test_blade.var --verbose 

loads: ## assign nodal forces to the joined mesh for use in FEA
	b3p_add_load $(YML) $(TMP)/joinedup_mesh.vtu

ccx: ## export a grid to ccx format, solve and postprocess (requires ccx to be installed and be in PATH)
	b3p_export_ccx $(TMP)/joinedup_mesh.vtu --out $(TMP)/ccxmodel.inp	
	ccx $(TMP)/ccxmodel_lc_forward_flap
	ccx $(TMP)/ccxmodel_lc_forward_edge
	b3p_frd2vtu $(TMP)/ccxmodel_lc_forward_edge.frd --out $(TMP)/ccx_fe.vtu
	b3p_frd2vtu $(TMP)/ccxmodel_lc_forward_flap.frd --out $(TMP)/ccx_ff.vtu


anbaprep: ## translate the section meshes to xdmf format
	python3 ../b3p/anba4_prep.py $(TMP)/msec_*.vtp


dockprep: ## prepare anba4 docker container
	docker build -t anba -f anba_docker .

andocker: ## solve using anba4 in a docker container
	docker run -v $(PWD)/..:/work -w /work anba "python3 b3p/anba4_solve.py examples/$(TMP)/msec_*.xdmf examples/$(TMP)/material_map.json --debug"


anbasolve: ## solve using anba4
	python3 ../b3p/anba4_solve.py $(TMP)/msec_*.xdmf $(TMP)/material_map.json --debug

all: geom mesh plybook drape combine mesh2d loads ccx ## run all targets in sequence
