
# b3p clean blade_test.yml
# tree temp_blade_portable 
# b3p build blade_test.yml
# tree temp_blade_portable
# b3p 

poetry run --project ~/projects/b3p b3p clean blade_test.yml 
poetry run --project ~/projects/b3p b3p build blade_test.yml geometry
poetry run --project ~/projects/b3p b3p build blade_test.yml mesh
poetry run --project ~/projects/b3p b3p build blade_test.yml drape
poetry run --project ~/projects/b3p b3p build blade_test.yml mass
poetry run --project ~/projects/b3p b3p build blade_test.yml apply-loads

poetry run --project ~/projects/b3p b3p ccx blade_test.yml 

poetry run --project ~/projects/b3p b3p 2d blade_test.yml 
