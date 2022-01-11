

geom:
	build_blade_geometry blade_test.yml

mesh:
	build_blade_structure blade_test.yml

plybook:
	build_plybook blade_test.yml --out temp/plybook.pck

drape:
	drape_mesh temp/plybook.pck temp/test_blade_web.vtp --key shell  --out temp/test_blade_shell.vtu 
	drape_mesh temp/plybook.pck temp/test_blade_w0.vtp --key web0  --out temp/test_blade_w0.vtu 

combine:
	combine_meshes temp/test_blade_shell.vtu  temp/test_blade_w0.vtu --out temp/joinedup_mesh.vtu

mesh2d:
	mesh_2d --vtu temp/joinedup_mesh.vtu --r "[5,20,35,50,65,80,90]" --var temp/test_blade.var


all: geom mesh plybook drape combine mesh2d
