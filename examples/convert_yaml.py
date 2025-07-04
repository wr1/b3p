#!/usr/bin/env python3

import yaml
import argparse
import sys


def convert_airfoil_yaml(input_file, output_file):
    """Load YAML and rewrite it using standard PyYAML formatting."""
    # Load input YAML
    try:
        with open(input_file, "r") as f:
            data = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: Input file {input_file} not found")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error: Invalid YAML in {input_file}: {e}")
        sys.exit(1)

    # Save output YAML with standard formatting
    try:
        with open(output_file, "w") as f:
            yaml.safe_dump(data, f, default_flow_style=None)
        print(f"Successfully wrote converted YAML to {output_file}")
    except Exception as e:
        print(f"Error: Failed to write output file {output_file}: {e}")
        sys.exit(1)


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Rewrite YAML file using standard PyYAML formatting."
    )
    parser.add_argument(
        "-i",
        "--input-yaml",
        type=str,
        default="S98_expanded.yml",
        help="Path to the input YAML file (default: S98_expanded.yml)",
    )
    parser.add_argument(
        "-o",
        "--output-yaml",
        type=str,
        default="converted_airfoil.yml",
        help="Path to the output YAML file (default: converted_airfoil.yml)",
    )
    args = parser.parse_args()

    convert_airfoil_yaml(args.input_yaml, args.output_yaml)


if __name__ == "__main__":
    main()
