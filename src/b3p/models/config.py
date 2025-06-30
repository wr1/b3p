"""Pydantic models for blade configuration."""

import numpy as np  # For vectorization of arrays
from pydantic import BaseModel, Field, validator, root_validator
from pathlib import Path
from typing import Dict, List, Optional, Union, Any

class IsotropicMaterial(BaseModel):
    """Model for isotropic materials."""

    E: float  # Young's modulus
    nu: float  # Poisson's ratio
    rho: float  # Density
    G: Optional[float] = None  # Shear modulus (optional)
    name: str = ""  # Material name


class AnisotropicMaterial(BaseModel):
    """Model for anisotropic materials."""

    Ex: float  # Young's modulus in x-direction
    Ey: Optional[float] = None  # Young's modulus in y-direction
    Ez: Optional[float] = None  # Young's modulus in z-direction
    Gxy: Optional[float] = None  # Shear modulus xy
    Gxz: Optional[float] = None  # Shear modulus xz
    Gyz: Optional[float] = None  # Shear modulus yz
    nu12: Optional[float] = None  # Poisson's ratio xy
    nu13: Optional[float] = None  # Poisson's ratio xz
    nu23: Optional[float] = None  # Poisson's ratio yz
    rho: float  # Density
    name: str = ""  # Material name


class PuckMaterial(AnisotropicMaterial):
    # Suggested Puck properties for failure criteria (add as needed):
    Xt: Optional[float] = None  # Tensile strength in x-direction
    Xc: Optional[float] = None  # Compressive strength in x-direction
    Yt: Optional[float] = None  # Tensile strength in y-direction
    Yc: Optional[float] = None  # Compressive strength in y-direction
    Zt: Optional[float] = None  # Tensile strength in z-direction
    Zc: Optional[float] = None  # Compressive strength in z-direction
    S12: Optional[float] = None  # Shear strength in 12 plane
    S13: Optional[float] = None  # Shear strength in 13 plane
    S23: Optional[float] = None  # Shear strength in 23 plane
    THETAF: Optional[float] = None  # Friction angle or similar
    MGF: Optional[float] = None  # Material growth factor
    ANU12: Optional[float] = None  # Additional Poisson's ratio
    ANU12f: Optional[float] = None  # Fiber-related Poisson's ratio
    E11_puck: Optional[float] = None  # Puck-specific modulus
    E11f: Optional[float] = None  # Fiber modulus for Puck


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
    splitstack: Optional[List[float]] = None  # e.g., [0.5, 0.5] if present


class LaminateConfig(BaseModel):
    """Laminate configuration."""

    slabs: Dict[str, Slab] = Field(default_factory=dict, description="Slab definitions")


class MaterialConfig(BaseModel):
    """Material configuration."""

    path: Optional[str] = None
    materials: Optional[Dict[str, Dict]] = None

    @root_validator(skip_on_failure=True)
    def check_materials(cls, values):
        """Ensure either path or materials is provided; treat input dict as materials if needed."""
        if 'path' not in values and 'materials' not in values:
            if len(values) > 0:
                values['materials'] = values  # Assign input dict to materials
            else:
                raise ValueError("Either materials.path or materials.materials must be provided")
        return values


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
        return v  # Root validator in MaterialConfig handles the check
