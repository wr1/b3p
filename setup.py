import setuptools

setuptools.setup(
    name="b3p",
    packages=["b3p", "utils"],
    install_requires=["fire"],
    entry_points={
        "console_scripts": [
            "b3p=b3p.b3p_cli:main",
            "b3p_slab_plot=utils.slab_plot:main",
            "b3p_drape_plot=utils.drape_plot:main",
            "b3p_aero=b3p.ccblade_run:main",
        ],
    },
)
