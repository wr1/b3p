# Usage Guide

## Basic Usage

B3P provides a command-line interface for running analyses:

```bash
b3p ccblade your_blade.yml
```

## Blade Definition

Blades are defined using YAML files. Here's a basic example:

```yaml
general:
  prefix: my_blade
  workdir: ./output

aero:
  bem:
    rated_power: 10e6  # 10 MW
    B: 3  # Number of blades
    rho: 1.225  # Air density
    polars:
      0.18: path/to/polar_18.dat
      0.21: path/to/polar_21.dat
      0.25: path/to/polar_25.dat
    max_tipspeed: 95.0
    uinf: [3, 5, 7, 9, 10, 11, 12, 13, 16, 20]
```

## Running CCBlade Analysis

The CCBlade analysis can be run programmatically:

```python
from b3p.ccblade_run import ccblade_run

# Run the analysis
analysis = ccblade_run("your_blade.yml")
```

## Interpreting Results

After running the analysis, results are saved to:

- `ccblade_output.csv`: Performance data at different wind speeds
- `ccblade_output.png`: Plot of performance metrics
- `polars.png`: Plot of interpolated airfoil polars
- `bladeloads.png`: Plot of blade loads

## Advanced Usage

### Custom Control Optimization

You can customize the control optimization process:

```python
from b3p.ccblade_run import ccblade_run, controloptimize

# Load the blade
blade = ccblade_run("your_blade.yml")

# Get the rotor
rotor = blade._setup_rotor()

# Custom control optimization
copt = controloptimize(
    rotor,
    max_tipspeed=95.0,
    rtip=blade.rtip,
    rating=15e6,  # Custom rating
    uinf=np.array([3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]),
    workdir="./custom_output"
)

# Optimize with custom parameters
copt.control_opt_below_rated(starting_tsr=8, starting_pitch=2)
output = copt.control_opt_above_rated()
``` 