"""Pydantic models for blade configuration."""

import numpy as np
from pydantic import (
    BaseModel,
    Field,
    field_validator,
    model_validator,
    # validator,
    # root_validator,
)
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
from .materials import IsotropicMaterial, AnisotropicMaterial, PuckMaterial


class GeneralConfig(BaseModel):
    """General configuration settings."""

    workdir: str = Field(default="output", description="Working directory path")
    prefix: str = Field(default="b3p", description="File prefix for outputs")


class Airfoil(BaseModel):
    """Airfoil data, either as a path or embedded coordinates."""

    path: Optional[str] = None
    name: Optional[str] = None
    xy: Optional[List[List[float]]] = None

    @field_validator("xy", mode="before")
    def validate_xy_path(cls, v):
        if isinstance(v, str):
            # If a path is provided, read the coordinates from the file
            file_path = Path(v)
            if not file_path.exists():
                raise ValueError(f"Airfoil file {v} does not exist")
            with open(file_path, "r") as f:
                lines = f.readlines()
            v = [list(map(float, line.split())) for line in lines if line.strip()]
        return v

    # @validator("xy", pre=True)
    # def validate_xy(cls, v):
    #     if v is not None:
    #         for point in v:
    #             if not isinstance(point, list) or len(point) != 2:
    #                 raise ValueError("xy must be a list of [x, y] coordinates")
    #     return v


class BemConfig(BaseModel):
    polars: Dict[float, str] = Field(
        default_factory=dict,
        description="Dictionary of polars with float keys and file paths as values",
    )
    rated_power: float = Field(
        default=10e6,
        description="Rated power of the turbine in watts",
    )
    B: int = Field(
        default=3,
        description="Number of blades",
    )
    rho: float = Field(
        default=1.225,
        description="Air density in kg/m^3",
    )
    tilt: float = Field(
        default=5.0,
        description="Tilt angle in degrees",
    )
    precone: float = Field(
        default=3.0,
        description="Precone angle in degrees",
    )
    shearExp: float = Field(
        default=0.1,
        description="Shear exponent for wind profile",
    )
    hubHt: float = Field(
        default=140.0,
        description="Hub height in meters",
    )
    mu: float = Field(
        default=1.81206e-5,
        description="Dynamic viscosity of air in kg/(m·s)",
    )
    yaw: float = Field(
        default=0.0,
        description="Yaw angle in degrees",
    )
    max_tipspeed: float = Field(
        default=95.0,
        description="Maximum tip speed in m/s",
    )
    uinf: List[float] = Field(
        default_factory=lambda: [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 20],
        description="List of inflow velocities in m/s",
    )


class AeroConfig(BaseModel):
    """Aerodynamic configuration."""

    airfoils: Dict[float, Airfoil] = Field(default_factory=dict)  # Keys as floats
    bem: BemConfig = Field(
        default_factory=BemConfig,
        description="BEM solver configuration",
    )


class MeshConfig(BaseModel):
    """Mesh configuration."""

    radii: str
    webs: Dict[str, Dict]
    panel_mesh_scale: List[List[int]] = Field(
        default_factory=lambda: [[0, 1]],
        description="Scale for panel mesh, e.g., [[0, 1]]",
    )
    n_web_points: int = Field(
        default=10,
        description="Number of points along the web",
    )
    n_chordwise_points: int = Field(
        default=100,
        description="Number of chordwise points in the mesh",
    )
    coordinates: Dict[str, Dict[str, Union[str, List[List[float]]]]] = Field(
        default_factory=lambda: {
            "d_te_offset": {
                "base": "d_te",
                "points": [[0, -0.1], [0.30, -0.2], [0.7, -0.3]],
            }
        },
        description="Coordinate offsets for the mesh",
    )
    bondline: Dict[str, Union[str, float, List[List[float]]]] = (
        Field(  # Updated to handle lists
            default_factory=lambda: {"type": "default", "thickness": 0.01, "width": []},
            description="Bondline settings",
        )
    )


class Slab(BaseModel):
    """Definition of a slab in the laminate configuration."""

    material: str
    cover: Dict[str, List[float]]  # e.g., {"d_w0": [-0.5, 0.5, 0]}
    slab: List[List[float]]  # e.g., [[0.03, 0], [0.10, 58]]
    ply_thickness: float
    key: List[Union[int, int]]  # e.g., [100, 2000]
    increment: List[Union[int, int]] = [1, -1]  # e.g., [1, -1]
    grid: str = "shell"  # e.g., "shell"
    splitstack: Optional[List[float]] = [1, 0]  # e.g., [0.5, 0.5] if present
    chamfers: Optional[List[Dict[str, Any]]] = None  # e.g., [[0.1, 0.2], [0.3, 0.4]]
    draping: Optional[str] = "plies"  # e.g., "blocks" or "plies"


class LaminateConfig(BaseModel):
    """Laminate configuration."""

    slabs: Dict[str, Slab] = Field(default_factory=dict, description="Slab definitions")


class BladeConfig(BaseModel):
    """Top-level blade configuration."""

    general: GeneralConfig = Field(default_factory=GeneralConfig)
    aero: AeroConfig = Field(default_factory=AeroConfig)
    mesh: MeshConfig = Field(default_factory=MeshConfig)
    laminates: LaminateConfig = Field(default_factory=LaminateConfig)
    materials: Dict[
        str, Union[IsotropicMaterial, AnisotropicMaterial, PuckMaterial]
    ] = Field(default_factory=dict)
    loads: Optional[Dict[str, Any]] = None
    damage: Optional[Dict[str, Any]] = None
    mesh2d: Optional[Dict[str, Any]]

    class Config:
        extra = "allow"

    @model_validator(mode="before")
    def validate_materials(cls, values):
        if not values.get("materials"):
            raise ValueError("materials must be provided in the configuration")
        return values

        # materials = values.get("materials")
        # if materials:
        #     for key, material in materials.items():
        #         if isinstance(material, IsotropicMaterial):
        #             material.type = "isotropic"
        #         elif isinstance(material, AnisotropicMaterial):
        #             material.type = "anisotropic"
        #         elif isinstance(material, PuckMaterial):
        #             material.type = "puck"
        #         else:
        #             raise ValueError(f"Unknown material type for {key}")
        # return values

    # @root_validator(skip_on_failure=True)
    # def check_materials(cls, values):
    #     if not values.get("materials"):
    #         raise ValueError("materials must be provided")
    #     return values
