import setuptools

setuptools.setup(
    name="b3p",
    packages=["b3p"],
    entry_points={
        "console_scripts": [
            "build_blade_geometry=b3p.build_blade_geometry:main",
            "build_blade_structure=b3p.build_blade_structure:main",
            "build_plybook=b3p.build_plybook:main",
            "drape_mesh=b3p.drape_mesh:main",
            "combine_meshes=b3p.combine_meshes:main",
            "mesh_2d=b3p.mesh_2d:main",
        ],
    },
)
