# Contributing

Thank you for your interest in contributing to B3P! This document provides guidelines and instructions for contributing.

## Development Setup

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/yourusername/b3p.git
   cd b3p
   ```
3. Install development dependencies:
   ```bash
   pip install -e ".[dev]"
   ```

## Code Style

We follow PEP 8 and use Black for code formatting:

```bash
black b3p
```

## Documentation

We use NumPy-style docstrings for all functions and classes. For example:

```python
def function(param1, param2):
    """
    Short description of the function.
    
    Parameters
    ----------
    param1 : type
        Description of param1
    param2 : type
        Description of param2
        
    Returns
    -------
    type
        Description of return value
    """
    # Function implementation
```

To build the documentation locally:

```bash
mkdocs serve
```

Then visit http://localhost:8000 to view the documentation.

## Testing

We use pytest for testing:

```bash
pytest
```

## Pull Request Process

1. Create a new branch for your feature or bugfix
2. Make your changes
3. Add tests for your changes
4. Update documentation if necessary
5. Run the tests to ensure they pass
6. Submit a pull request

## Code of Conduct

Please be respectful and inclusive in your interactions with others in the community. 