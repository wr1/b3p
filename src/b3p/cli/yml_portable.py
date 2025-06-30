import logging
from pathlib import Path
from ruamel.yaml import YAML
from importlib.resources import files
from ..models.config import BladeConfig, Airfoil
import numpy as np
import re
import os

logger = logging.getLogger(__name__)

def load_airfoil(af: str, prefix: str = "") -> tuple[Optional[str], Optional[list]]:
    """Load an airfoil from an xfoil formatted text file."""
    af_path = os.path.join(prefix, af) if prefix else af
    if not os.path.exists(af_path):
        logger.warning(f"Airfoil file {af_path} not found, keeping as reference")
        return None, None
    with open(af_path, "r") as f:
        name = f.readline().strip()
        floats = re.findall(r"\d+\.\d+", name)
        if len(floats) == 2:
            offset = 0
            name = Path(af).stem
        elif name == "":
            offset = 1
            name = Path(af).stem
        else:
            offset = 1
    return name, np.loadtxt(af_path, skiprows=offset).tolist()

def load_airfoils(afs: dict, prefix: str = "") -> dict:
    """Load a list of airfoils from a dictionary."""
    dct = {}
    for key, value in afs.items():
        if isinstance(value, dict) and "xy" in value:
            dct[key] = value
        else:
            name, xy = load_airfoil(value, prefix)
            if name is None:
                dct[key] = {"path": value}
            else:
                dct[key] = {"xy": xy, "name": name, "path": value}
                logger.info(f"Imported airfoil {name} at thickness {key}")
    return dct

def yaml_make_portable(yaml_file: Path) -> BladeConfig:
    """Load and validate a YAML file, merging with defaults."""
    yaml = YAML()
    prefix = str(yaml_file.parent)

    # Load user-provided YAML
    try:
        user_data = yaml.load(yaml_file.open("r"))
    except Exception as e:
        logger.error(f"Failed to load YAML {yaml_file}: {e}")
        raise

    # Load default YAML
    default_path = files("b3p.models") / "defaults.yaml"
    with default_path.open("r") as f:
        default_data = yaml.load(f)

    # Merge defaults with user data (user data overrides defaults)
    merged_data = default_data.copy()
    merged_data.update(user_data)

    # Process airfoils
    if "aero" in merged_data and "airfoils" in merged_data["aero"]:
        merged_data["aero"]["airfoils"] = load_airfoils(merged_data["aero"]["airfoils"], prefix)

    # Process subsection files (materials, loads, laminates)
    subsections = ["materials", "loads", "laminates"]
    for s in subsections:
        if s in merged_data and isinstance(merged_data[s], str) and merged_data[s].endswith(".yml"):
            sub_path = os.path.join(prefix, merged_data[s])
            if os.path.exists(sub_path):
                logger.info(f"Loading {s} from: {sub_path}")
                merged_data[s] = yaml.load(open(sub_path, "r"))
            else:
                logger.warning(f"File {sub_path} not found, keeping as reference")
                merged_data[s] = {"path": merged_data[s]}

    # Validate with Pydantic
    try:
        config = BladeConfig(**merged_data)
    except Exception as e:
        logger.error(f"YAML validation failed: {e}")
        raise

    logger.info(f"Loaded and validated YAML from {yaml_file}")
    return config

def save_yaml(of: str, config: BladeConfig):
    """Save BladeConfig to YAML file."""
    yaml = YAML()
    yaml.default_flow_style = True
    with open(of, "w") as f:
        yaml.dump(config.dict(exclude_none=True), f)
        logger.info(f"Written to: {of}")
