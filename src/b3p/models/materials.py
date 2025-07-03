"""Pydantic models for materials."""

from pydantic import BaseModel
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
