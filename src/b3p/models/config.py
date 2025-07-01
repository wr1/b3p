"""Pydantic models for blade configuration."""

import numpy as np
from pydantic import BaseModel, Field, validator, root_validator
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

    @validator("xy", pre=True)
    def validate_xy(cls, v):
        if v is not None:
            for point in v:
                if not isinstance(point, list) or len(point) != 2:
                    raise ValueError("xy must be a list of [x, y] coordinates")
        return v


class AeroConfig(BaseModel):
    """Aerodynamic configuration."""

    airfoils: Dict[float, Airfoil] = Field(default_factory=dict)  # Keys as floats


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
    cover: Dict[str, List[Union[float, int]]]  # e.g., {"d_w0": [-0.5, 0.5, 0]}
    slab: List[List[Union[float, int]]]  # e.g., [[0.03, 0], [0.10, 58]]
    ply_thickness: float
    key: List[Union[int, float]]  # e.g., [100, 2000]
    increment: List[Union[int, int]] = [1, -1]  # e.g., [1, -1]
    grid: str  # e.g., "shell"
    splitstack: Optional[List[float]] = [1, 0]  # e.g., [0.5, 0.5] if present
    chamfers: Optional[List[Dict[str, Any]]] = None  # e.g., [[0.1, 0.2], [0.3, 0.4]]


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

    @root_validator(skip_on_failure=True)
    def check_materials(cls, values):
        if not values.get("materials"):
            raise ValueError("materials must be provided")
        return values
