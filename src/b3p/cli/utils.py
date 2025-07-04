import subprocess


def wslpath_convert(path: str, to_windows: bool = False) -> str:
    """Converts paths between WSL and Windows formats using wslpath."""
    command = ["wslpath", "-w" if to_windows else "-u", path]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    return result.stdout.strip()


def check_existing_outputs(files):
    """Check if all specified files exist."""
    return all(f.exists() for f in files)
