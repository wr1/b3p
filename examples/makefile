
# https://betterprogramming.pub/how-to-create-a-self-documenting-makefile-533ebf8f82e2
.PHONY: help
help: ## Show this help
	@egrep -h '\s##\s' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'


YML = blade_test_portable.yml
TMP= temp_b50_portable

portable: ## make self contained yaml file 
	b3p_yml_portable blade_test.yml 

clean: ## remove all temporary files
	rm -rf $(TMP)

geom: ## create blade shell geometry (no shearwebs)
	b3p_blade_geometry $(YML) 

mesh: ## create blade structure 
	b3p_blade_structure $(YML)

plybook: ## create a plybook out of slab definition
	b3p_plybook $(YML) --out $(TMP)/plybook.pck

drape: ## drape the shell and the shear web
	b3p_drape $(TMP)/plybook.pck $(TMP)/test_blade_web.vtp --key shell  --out $(TMP)/test_blade_web.vtu 
	b3p_drape $(TMP)/plybook.pck $(TMP)/test_blade_w4.vtp --key web4  --out $(TMP)/test_blade_w04.vtu 
	b3p_drape $(TMP)/plybook.pck $(TMP)/test_blade_w3.vtp --key web3  --out $(TMP)/test_blade_w03.vtu 
	b3p_drape $(TMP)/plybook.pck $(TMP)/test_blade_w5.vtp --key web5  --out $(TMP)/test_blade_w05.vtu 

combine: ## add the meshes together 
	b3p_meshcombine $(TMP)/test_blade_web.vtu $(TMP)/test_blade_w05.vtu $(TMP)/test_blade_w04.vtu $(TMP)/test_blade_w03.vtu --out $(TMP)/joinedup_mesh.vtu

# $(TMP)/test_blade_w04.vtu
mesh2d: ## extract 2D meshes 
	b3p_2dmesh --vtu $(TMP)/joinedup_mesh.vtu --r "[2,10,20,30,40,50,60,70,80,90]" --var $(TMP)/test_blade.var --verbose 

loads: ## assign nodal forces to the joined mesh for use in FEA
	b3p_add_load $(YML) $(TMP)/joinedup_mesh.vtu

ccxpre: ## export a grid to ccx format
	b3p_export_ccx $(TMP)/joinedup_mesh.vtu --out $(TMP)/ccxmodel.inp	--matmap $(TMP)/material_map.json --merge_adjacent --quadratic --force_iso

ccxsolve: ## run calculix on the ccx output models, and postprocess (requires ccx to be installed and be in PATH)
	ccx $(TMP)/ccxmodel_lc_forward_flap
	ccx $(TMP)/ccxmodel_lc_forward_edge

ccxpost: ## postprocess the ccx results
	b3p_frd2vtu $(TMP)/ccxmodel_lc_forward_edge.frd --out $(TMP)/ccx_fe.vtu
	b3p_frd2vtu $(TMP)/ccxmodel_lc_forward_flap.frd --out $(TMP)/ccx_ff.vtu

ccx: ccxpre ccxsolve ccxpost ## run ccxpre and ccxsolve

anbaprep: ## translate the section meshes to xdmf format
	b3p_anbaprep $(TMP)/msec_*.vtp

dockerprep: ## prepare anba4 docker container
	docker build -t anba -f anba_docker .

anbasolve: ## solve using anba4 in a docker container
	docker run -v $(PWD)/..:/work -w /work anba "python3 b3p/anba4_solve.py examples/$(TMP)/msec_*.xdmf examples/$(TMP)/material_map.json --debug"

all: geom mesh plybook drape combine mesh2d loads ccx anbaprep anbasolve ## run all targets in sequence
