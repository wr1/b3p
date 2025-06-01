import subprocess
from pathlib import Path

def wslpath_convert(path: str, to_windows: bool = False) -> str:
    """
    Converts paths between WSL and Windows formats using wslpath.

    :param path: The path to convert.
    :param to_windows: Set to True to convert WSL to Windows; False for Windows to WSL.
    :return: The converted path as a string.
    """
    command = ["wslpath", "-w" if to_windows else "-u", path]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    return result.stdout.strip()

def check_existing_outputs(files):
    """
    Check if all specified files exist.

    :param files: List of Path objects to check.
    :return: True if all files exist, False otherwise.
    """
    return all(f.exists() for f in files)
