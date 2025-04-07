from .geometry import blade
import importlib.metadata

NAME = "b3p"

try:
    __version__ = importlib.metadata.version(NAME)
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"
