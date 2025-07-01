"""Material database class for loading, managing, and exporting material properties."""

import ruamel.yaml
from pathlib import Path
from pydantic import BaseModel
from typing import Dict, Any, Union, Optional
from ..materials.aniso_material import AnisotropicMaterial
from ..materials.iso_material import IsotropicMaterial


class MaterialDatabase:
    """Database for materials, supporting YAML loading and property export."""

    def __init__(self) -> None:
        self.materials: Dict[str, Union[IsotropicMaterial, AnisotropicMaterial]] = {}

    def load_from_yaml(self, file_path: str) -> None:
        """Load materials from a YAML file."""
        path = Path(file_path)
        with path.open("r") as file:
            data = ruamel.yaml.YAML().load(file)
            for key, props in data.items():
                if "name" not in props:
                    props["name"] = key
                # Determine material type based on properties (simple heuristic)
                if "E" in props and "nu" in props and "rho" in props:
                    self.materials[key] = IsotropicMaterial(**props)  # type: ignore
                else:
                    self.materials[key] = AnisotropicMaterial(**props)  # type: ignore

    def write_to_yaml(self, file_path: str) -> None:
        """Write materials to a YAML file."""
        data: Dict[str, Dict[str, Any]] = {
            key: mat.dict() for key, mat in self.materials.items()
        }
        path = Path(file_path)
        with path.open("w") as file:
            ruamel.yaml.YAML().dump(data, file)

    def get_properties_ccx(self, material_key: str) -> Dict[str, Any]:
        """Supply properties mapped to CCX (CalculiX) naming scheme."""
        material = self.materials.get(material_key)
        if not material:
            return {}
        if isinstance(material, IsotropicMaterial):
            return {
                "YOUNG": material.E,
                "POISS": material.nu,
                "DENS": material.rho,
                "SHEAR": material.G,
            }
        else:  # Anisotropic
            return {
                "E1": material.Ex,
                "E2": material.Ey,
                "E3": material.Ez,
                "NU12": material.nu12,
                "NU13": material.nu13,
                "NU23": material.nu23,
                "G12": material.Gxy,
                "G13": material.Gxz,
                "G23": material.Gyz,
                "DENS": material.rho,
            }

    def get_properties_anba(self, material_key: str) -> Dict[str, Any]:
        """Supply properties mapped to ANBA naming scheme."""
        material = self.materials.get(material_key)
        if not material:
            return {}
        if isinstance(material, IsotropicMaterial):
            return {
                "MODULUS": material.E,
                "POISSON": material.nu,
                "DENSITY": material.rho,
                "SHEARMOD": material.G,
            }
        else:  # Anisotropic
            return {
                "E_X": material.Ex,
                "E_Y": material.Ey,
                "E_Z": material.Ez,
                "POISSON_XY": material.nu12,
                "POISSON_XZ": material.nu13,
                "POISSON_YZ": material.nu23,
                "SHEAR_XY": material.Gxy,
                "SHEAR_XZ": material.Gxz,
                "SHEAR_YZ": material.Gyz,
                "DENSITY": material.rho,
            }

    def is_puck_material(self, material_key: str) -> bool:
        """Check if the material has a full set of Puck properties."""
        material = self.materials.get(material_key)
        if isinstance(material, AnisotropicMaterial):
            puck_props = [
                "Xt",
                "Xc",
                "Yt",
                "Yc",
                "Zt",
                "Zc",
                "S12",
                "S13",
                "S23",
                "THETAF",
                "MGF",
                "ANU12",
                "ANU12f",
                "E11_puck",
                "E11f",
            ]
            return all(getattr(material, prop, None) is not None for prop in puck_props)
        return False  # Only anisotropic materials can be Puck materials

    def add_material(
        self, key: str, material: Union[IsotropicMaterial, AnisotropicMaterial]
    ) -> None:
        """Add a new material to the database."""
        self.materials[key] = material


# Example usage (not executed here)
# db = MaterialDatabase()
# db.load_from_yaml('materials_si.yml')
# props = db.get_properties_ccx('balsa110')
# print(props)
# if db.is_puck_material('some_material'):
#     print("This is a Puck material")
