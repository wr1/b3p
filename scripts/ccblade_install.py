import subprocess
from pathlib import Path


def build_ccblade_dependency():
    repo_name = "ccblade"
    build_dir = Path("build") / repo_name
    if not build_dir.exists():
        subprocess.run(
            ["git", "clone", "https://github.com/WISDEM/CCBlade.git", build_dir],
            check=True,
        )
    subprocess.run(["pip", "install", "."], cwd=build_dir, check=True)


if __name__ == "__main__":
    build_ccblade_dependency()
