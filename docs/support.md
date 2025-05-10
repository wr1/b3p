# Support

If you encounter issues, have questions, or want to suggest features for **B3P** (Blade Preprocessor), this guide outlines how to seek help.

## GitHub Repository

The primary support channel is the B3P GitHub repository: [https://github.com/wr1/b3p](https://github.com/wr1/b3p).

### Filing Issues

For bugs, installation problems, or feature requests:

1. Check existing issues to avoid duplicates.
2. Open a new issue with:

   - A clear title (e.g., "Error in 2D mesh generation").
   - Description of the problem, including error messages and steps to reproduce.
   - Environment details (OS, Python version, B3P version, dependencies).
   - Relevant logs or screenshots.

### Discussions

Use GitHub Discussions for general questions, ideas, or community feedback:

- Ask about usage, best practices, or clarification on documentation.
- Share ideas for new features or improvements.

## Troubleshooting Common Issues

- **Installation failures**: Ensure Python 3.8+, `pip`, and dependencies are up-to-date. Use a virtual environment to avoid conflicts.
- **ANBA4 errors**: Verify the `anba4-env` Conda environment is activated and ANBA4 is installed.
- **CalculiX not found**: Confirm `ccx` is in your PATH or specify with `--ccxexe`.
- **YAML errors**: Check file paths and syntax. Use `b3p yml_portable` to embed linked files.
- **Memory issues**: Reduce radial positions (`mesh.radii`) or use a machine with more RAM.

## Community and Resources

- **Documentation**: Refer to [Usage](usage.md), [Input File Format](use/inputfile.md), and [Examples](examples/blade_test.md) for guidance.
- **Examples**: The `examples` folder in the repository includes workflows to learn from.
- **Paraview**: Use Paraview to visualize `.vtp`, `.vtu`, and `.xdmf` output files.


