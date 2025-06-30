from pydantic import BaseModel
from typing import Optional


class IsotropicMaterial(BaseModel):
    """Model for isotropic materials."""

    E: float  # Young's modulus
    nu: float  # Poisson's ratio
    rho: float  # Density
    G: Optional[float] = None  # Shear modulus (optional)
    name: str = ""  # Material name
