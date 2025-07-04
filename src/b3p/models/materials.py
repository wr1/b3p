"""Pydantic models for materials."""

from pydantic import BaseModel, Field
from typing import Optional

class IsotropicMaterial(BaseModel):
    """Model for isotropic materials."""

    E: float  # Young's modulus
    nu: float  # Poisson's ratio
    rho: float  # Density
    G: Optional[float] = None  # Shear modulus (optional)
    name: str = ""  # Material name

class AnisotropicMaterial(BaseModel):
    """Model for anisotropic materials."""

    # remove aliases
    Ex: float #= Field(aliases=["e11"])  # Young's modulus in x-direction
    Ey:   Optional[float] = Field(aliases=["e22"], default=None)  # Young's modulus in y-direction
    Ez:   Optional[float] = Field(aliases=["e33"], default=None)  # Young's modulus in z-direction
    Gxy:  Optional[float] = Field(aliases=["g12"], default=None)  # Shear modulus xy
    Gxz:  Optional[float] = Field(aliases=["g13"], default=None)  # Shear modulus xz
    Gyz:  Optional[float] = Field(aliases=["g23"], default=None)  # Shear modulus yz
    nuxy: Optional[float] = Field(aliases=["nu12"], default=None)  # Poisson's ratio xy
    nuxz: Optional[float] = Field(aliases=["nu13"], default=None)  # Poisson's ratio xz
    nuyz: Optional[float] = Field(aliases=["nu23"], default=None)  # Poisson's ratio yz
    rho: float  # Density
    name: str = ""  # Material name

class PuckMaterial(AnisotropicMaterial):
    # Suggested Puck properties for failure criteria
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
