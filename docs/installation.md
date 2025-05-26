# Installation

This guide covers the installation of **B3P** (Blade Preprocessor) and its dependencies. B3P is tested on Ubuntu 24.04 and Windows WSL2 (do not use WSL1). A Python environment (version 3.8 or higher) is required.

## Prerequisites

Before installing B3P, ensure the following are installed:

- **Python 3.8+**: Install via your package manager (e.g., `apt install python3` on Ubuntu) or from [python.org](https://www.python.org).
- **pip**: Python package manager, typically included with Python.
- **Git**: For cloning the repository (e.g., `apt install git` on Ubuntu).
- **Conda** (optional): For managing dependencies like ANBA4 (install Miniconda or Anaconda from [conda.io](https://conda.io)).
- **CalculiX** (optional): Required for 3D FEA (install via `apt install calculix-ccx` on Ubuntu or follow [CalculiX instructions](http://www.calculix.de)).

## Installation Options

B3P can be installed either from PyPI or from source.

### Option 1: Install from PyPI

The simplest way to install B3P is using pip:

```bash
pip install b3p
```

This installs B3P and its core dependencies (e.g., `pyvista`, `numpy`, `pandas`).

### Option 2: Install from Source

To install from source (e.g., for development or access to the latest features):

1. Clone the repository:
   ```bash
   git clone https://github.com/wr1/b3p.git
   cd b3p
   ```

2. Install B3P in editable mode:
   ```bash
   pip install -e .
   ```

3. (Optional) Install development dependencies:
   ```bash
   pip install -e ".[dev]"
   ```

### Additional Dependencies

For specific analyses, additional tools may be required:

- **ANBA4**: For 2D sectional analysis, install ANBA4 in a Conda environment (e.g., `anba4-env`). Follow ANBA4's installation instructions or use the `default_install.sh` script in the B3P repository.
- **CCBlade**: For aerodynamic analysis, included with B3P but requires airfoil polar data.

If using the `default_install.sh` script:

```bash
cd b3p
bash default_install.sh
```

This script sets up ANBA4, CalculiX, and other dependencies in a Conda environment.

## Verification

Verify the installation by checking the B3P version:

```bash
b3p --version
```

Run the test suite to ensure functionality:

```bash
pytest
```

## Troubleshooting

- **pip errors**: Ensure `pip` is up-to-date (`pip install --upgrade pip`).
- **Conda issues**: Verify Conda is activated (`conda activate anba4-env`) and the environment includes ANBA4.
- **CalculiX not found**: Confirm `ccx` is in your PATH (`which ccx`) or specify the executable path with `--ccxexe`.
- **Permission errors**: Use `sudo` for system-wide installations or a virtual environment.

For further assistance, check the [GitHub issues page](https://github.com/wr1/b3p/issues).

