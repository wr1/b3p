# Contributing to B3P

Thank you for your interest in contributing to **B3P** (Blade Preprocessor)! This guide outlines the process for contributing code, documentation, or other improvements.

## Getting Started

### 1. Fork and Clone the Repository

1. Fork the B3P repository on GitHub: [https://github.com/wr1/b3p](https://github.com/wr1/b3p).
2. Clone your fork locally:
   ```bash
   git clone https://github.com/yourusername/b3p.git
   cd b3p
   ```

### 2. Set Up the Development Environment

1. Install B3P in editable mode with development dependencies:
   ```bash
   pip install -e ".[dev]"
   ```
2. (Optional) Set up a Conda environment for ANBA4:
   ```bash
   bash scripts/anba4_wsl_setup.sh
   ```

### 3. Create a Branch

Create a new branch for your feature or bug fix:

```bash
git checkout -b feature/your-feature-name
```

## Development Guidelines

### Code Style

- Follow **PEP 8** for Python code.
- Use **Black** for automatic code formatting:
   ```bash
   black b3p
   ```
- Use **NumPy-style docstrings** for functions and classes:
   ```python
   def example_function(param1: int, param2: str) -> bool:
       """
       Short description of the function.
   
       Parameters
       ----------
       param1 : int
           Description of param1.
       param2 : str
           Description of param2.
   
       Returns
       -------
       bool
           Description of return value.
       """
       pass
   ```

### Testing

- Write tests for new features or bug fixes using **pytest**.
- Run the test suite to ensure all tests pass:
   ```bash
   pytest
   ```
- If a test is commented out, assume it’s failing. Fix the underlying issue rather than re-enabling it.

### Documentation

- Update documentation in the `docs/` folder for new features or changes.
- Use Markdown for documentation files.
- Build and preview documentation locally:
   ```bash
   mkdocs serve
   ```
   Visit `http://localhost:8000` to view the documentation.

## Submitting Changes

### 1. Commit Changes

Write clear, concise commit messages:

```bash
git commit -m "Add feature X to module Y"
```

### 2. Push to Your Fork

Push your branch to your GitHub fork:

```bash
git push origin feature/your-feature-name
```

### 3. Create a Pull Request

1. Open a pull request (PR) from your branch to the `main` branch of the main repository.
2. Include a detailed description of your changes, including:
   - Purpose of the change.
   - Related issues (e.g., `Fixes #123`).
   - Any testing performed.
3. Ensure all tests pass and documentation is updated.

## Code of Conduct

Be respectful and inclusive in all interactions. Follow the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org).

## Troubleshooting

- **Test failures**: Check error messages and verify dependency versions.
- **Formatting issues**: Run `black` before committing.
- **Documentation errors**: Validate Markdown syntax and test with `mkdocs serve`.

For questions, reach out via GitHub issues or discussions.

