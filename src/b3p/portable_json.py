import json
from pathlib import Path
from abc import ABC, abstractmethod
from ruamel.yaml import YAML
import numpy as np
import logging
from datetime import datetime

logger = logging.getLogger(__name__)

class FileHandler(ABC):
    """Base class for handling different file types."""
    @abstractmethod
    def validate(self, file_path: Path) -> bool:
        """Validate the file format and contents."""
        pass

    @abstractmethod
    def load(self, file_path: Path) -> dict:
        """Load the file into a serializable dictionary."""
        pass

class YamlHandler(FileHandler):
    """Handler for YAML files."""
    def validate(self, file_path: Path) -> bool:
        if not file_path.exists():
            logger.error(f"YAML file {file_path} does not exist")
            return False
        if file_path.suffix not in {".yml", ".yaml"}:
            logger.error(f"File {file_path} is not a YAML file")
            return False
        return True

    def load(self, file_path: Path) -> dict:
        yaml = YAML()
        with open(file_path, "r") as f:
            return yaml.load(f)

class AirfoilHandler(FileHandler):
    """Handler for XFOIL-format airfoil files."""
    def validate(self, file_path: Path) -> bool:
        if not file_path.exists():
            logger.error(f"Airfoil file {file_path} does not exist")
            return False
        try:
            with open(file_path, "r") as f:
                first_line = f.readline().strip()
                floats = [float(x) for x in first_line.split() if x]
                return len(floats) == 2 or first_line == ""
        except Exception as e:
            logger.error(f"Invalid airfoil format in {file_path}: {e}")
            return False

    def load(self, file_path: Path) -> dict:
        with open(file_path, "r") as f:
            name = f.readline().strip()
            floats = [float(x) for x in name.split() if x]
            offset = 0 if len(floats) == 2 else 1
            if offset == 1:
                name = file_path.stem
        xy = np.loadtxt(file_path, skiprows=offset).tolist()
        return {"name": name, "xy": xy}

class PortableJsonConfig:
    """Manages creation and loading of portable JSON configurations."""
    def __init__(self):
        self.handlers = {
            "yaml": YamlHandler(),
            "airfoil": AirfoilHandler(),
        }
        self.dependency_keys = {
            "aero.airfoils": ("airfoil", "dict"),
            "materials": ("yaml", "file"),
            "loads": ("yaml", "file"),
            "laminates": ("yaml", "file"),
        }

    def register_handler(self, type_name: str, handler: FileHandler) -> None:
        """Register a custom file handler."""
        self.handlers[type_name] = handler

    def make_portable(self, yaml_file: Path, output_json: Path = None) -> Path:
        """Create a portable JSON from a YAML file and its dependencies."""
        yaml_file = Path(yaml_file)
        if output_json is None:
            output_json = yaml_file.parent / f"{yaml_file.stem}_portable.json"

        # Load and validate main YAML
        if not self.handlers["yaml"].validate(yaml_file):
            raise ValueError(f"Invalid YAML file: {yaml_file}")
        config = self.handlers["yaml"].load(yaml_file)

        # Initialize JSON structure
        portable_data = {
            "config": config,
            "dependencies": {},
            "metadata": {
                "source_yaml": str(yaml_file.name),
                "created": datetime.utcnow().isoformat() + "Z",
                "version": "1.0"
            }
        }

        # Process dependencies
        for key, (file_type, structure) in self.dependency_keys.items():
            if key not in config:
                continue
            handler = self.handlers.get(file_type)
            if not handler:
                logger.warning(f"No handler for file type {file_type} at {key}")
                continue

            portable_data["dependencies"][key] = {}
            if structure == "dict":
                # Handle dictionary of files (e.g., aero.airfoils)
                for af_key, af_path in config[key].items():
                    af_path = Path(yaml_file.parent) / af_path
                    if not handler.validate(af_path):
                        raise ValueError(f"Invalid file {af_path} for {key}")
                    portable_data["dependencies"][key][af_key] = handler.load(af_path)
            elif structure == "file":
                # Handle single file (e.g., materials)
                file_path = Path(yaml_file.parent) / config[key]
                if not handler.validate(file_path):
                    raise ValueError(f"Invalid file {file_path} for {key}")
                portable_data["dependencies"][key] = handler.load(file_path)

        # Update workdir to output directory
        output_dir = output_json.parent
        portable_data["config"]["general"]["workdir"] = str(output_dir)

        # Save JSON
        with open(output_json, "w") as f:
            json.dump(portable_data, f, indent=2)
        logger.info(f"Saved portable JSON to {output_json}")
        return output_json

    def load_portable(self, json_file: Path) -> dict:
        """Load a portable JSON and return the configuration dictionary."""
        json_file = Path(json_file)
        if not json_file.exists():
            raise ValueError(f"JSON file {json_file} does not exist")
        with open(json_file, "r") as f:
            portable_data = json.load(f)

        # Validate structure
        if "config" not in portable_data or "dependencies" not in portable_data:
            raise ValueError(f"Invalid portable JSON format in {json_file}")

        # Update workdir to current directory
        portable_data["config"]["general"]["workdir"] = str(json_file.parent)
        return portable_data["config"]

    def unpack(self, json_file: Path, output_dir: Path) -> Path:
        """Unpack a portable JSON into a YAML and separate files."""
        json_file = Path(json_file)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Load JSON
        with open(json_file, "r") as f:
            portable_data = json.load(f)

        # Save main YAML
        yaml = YAML()
        yaml.default_flow_style = True
        config = portable_data["config"]
        output_yaml = output_dir / "config.yml"
        with open(output_yaml, "w") as f:
            yaml.dump(config, f)
        logger.info(f"Saved YAML to {output_yaml}")

        # Save dependencies as files
        for key, (file_type, structure) in self.dependency_keys.items():
            if key not in portable_data["dependencies"]:
                continue
            handler = self.handlers.get(file_type)
            if not handler:
                continue

            if structure == "dict":
                for af_key, data in portable_data["dependencies"][key].items():
                    output_path = output_dir / f"{af_key}.txt"
                    handler.save(data, output_path)
                    config[key][af_key] = str(output_path.relative_to(output_dir))
            elif structure == "file":
                output_path = output_dir / f"{key.replace('.', '_')}.yml"
                handler.save(portable_data["dependencies"][key], output_path)
                config[key] = str(output_path.relative_to(output_dir))

        # Update workdir and save YAML again
        config["general"]["workdir"] = str(output_dir)
        with open(output_yaml, "w") as f:
            yaml.dump(config, f)

        logger.info(f"Unpacked portable JSON to {output_dir}")
        return output_dir

def make_portable(yaml_file: str) -> None:
    """Convenience function to create a portable JSON."""
    PortableJsonConfig().make_portable(Path(yaml_file))

def load_portable(json_file: str) -> dict:
    """Convenience function to load a portable JSON."""
    return PortableJsonConfig().load_portable(Path(json_file))
