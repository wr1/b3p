# https://betterprogramming.pub/how-to-create-a-self-documenting-makefile-533ebf8f82e2
.PHONY: help
help: ## Show this help
	@egrep -h '\s##\s' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

YML = blade_test.yml
TMP= temp_blade_portable

build: ## make self contained yaml file 
	b3p --yml=$(YML) clean build 

mass:
	b3p --yml=$(YML) mass
	
ccx: ccxprep ccxsolve ccxpost ## prepare and run ccx models

ccxprep: ## prepare ccx model
	b3p --yml=$(YML) ccxprep --buckling=False

ccxprepb: ## prepare ccx model for buckling run
	b3p --yml=$(YML) ccxprep --buckling=True

ccxsolve: ## run ccx on the model
	b3p --yml=$(YML) ccxsolve forward_flap 

ccxpost: ## postprocess ccx output
	b3p --yml=$(YML) ccxpost
	b3p --yml=$(YML) ccxplot

anba: ## run anba
	b3p --yml=$(YML) mesh2d  90
	
# --z_start=0.0001 -z_end=99.5 --nsec=50

aero: ## run aero analysis
	b3p_aero --blade=$(YML)

all: build ccx anba aero ## run all steps

clean: ## remove all temporary files
	rm -rf $(TMP)

anbasolve: # solve using anba4 in a docker container
	docker run -v $(PWD)/..:/work -w /work anba "python3 b3p/anba4_solve.py examples/$(TMP)/msec_*.xdmf examples/$(TMP)/material_map.json --debug"

all: build ccx anba anbasolve aero ## run all targets in sequence
