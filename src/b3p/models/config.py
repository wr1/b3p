"""Pydantic models for blade configuration."""

from pydantic import BaseModel, Field, validator
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
from ..materials.iso_material import IsotropicMaterial
from ..materials.aniso_material import AnisotropicMaterial

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

    bondline: Dict[str, Union[str, float, List[List[float]]]] = Field(  # Updated to handle lists
        default_factory=lambda: {"type": "default", "thickness": 0.01, "width": []},
        description="Bondline settings",
    )

class Slab(BaseModel):
    """Definition of a slab in the laminate configuration."""

    material: str
    cover: Dict[str, List[Union[float, int]]]  # e.g., {"d_w0": [-0.5, 0.5, 0]}
    slab: List[List[Union[float, int]]]  # e.g., [[0.03, 0], [0.10, 58]]
    ply_thickness: float
    key: List[Union[int, float]]  # e.g., [100, 2000]
    increment: List[Union[int, float]]  # e.g., [1, -1]
    grid: str  # e.g., "shell"
    splitstack: Optional[List[float]] = None  # e.g., [0.5, 0.5] if present

class LaminateConfig(BaseModel):
    """Laminate configuration."""

    slabs: Dict[str, Slab] = Field(default_factory=dict, description="Slab definitions")

class MaterialConfig(BaseModel):
    """Material configuration."""

    path: Optional[str] = None
    materials: Optional[Dict[str, Union[IsotropicMaterial, AnisotropicMaterial]]] = None

class BladeConfig(BaseModel):
    """Top-level blade configuration."""

    general: GeneralConfig = Field(default_factory=GeneralConfig)
    aero: AeroConfig = Field(default_factory=AeroConfig)
    mesh: MeshConfig = Field(default_factory=MeshConfig)
    laminates: LaminateConfig = Field(default_factory=LaminateConfig)
    materials: MaterialConfig = Field(default_factory=MaterialConfig)
    loads: Optional[Dict[str, Any]] = None
    damage: Optional[Dict[str, Any]] = None

    class Config:
        extra = "allow"

    @validator("materials")
    def validate_materials(cls, v):
        if v.path is None and v.materials is None:
            raise ValueError(
                "Either materials.path or materials.materials must be provided"
            )
        return v
