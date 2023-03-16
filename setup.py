import setuptools

setuptools.setup(
    name="b3p",
    packages=["b3p"],
    install_requires=["fire"],
    entry_points={
        "console_scripts": [
            "b3p=b3p.b3p_cli:main"
            # "b3p_blade_geometry=b3p.build_blade_geometry:main",
            # "b3p_blade_structure=b3p.build_blade_structure:main",
            # # "b3p_plybook=b3p.build_plybook:main",
            # "b3p_drape=b3p.drape_mesh:main",
            # # "b3p_meshcombine=b3p.combine_meshes:main",
            # "b3p_2dmesh=b3p.mesh_2d:main",
            # "b3p_add_load=b3p.add_load_to_mesh:main",
            # "b3p_export_ccx=b3p.mesh2ccx:main",
            # "b3p_frd2vtu=b3p.frd2vtu:main",
            # "b3p_anbaprep=b3p.anba4_prep:main",
            # "b3p_yml_portable=utils.yml_portable:main",
        ],
    },
)
