# Contributing to Phytovenomics

Thank you for considering contributing to the Phytovenomics project! This document provides guidelines and instructions for contributing to make the process smooth and effective for everyone involved.

## Table of Contents

1. [Code of Conduct](#code-of-conduct)
2. [Getting Started](#getting-started)
3. [Development Environment](#development-environment)
4. [Branching Strategy](#branching-strategy)
5. [Coding Standards](#coding-standards)
6. [Documentation](#documentation)
7. [Testing](#testing)
8. [Submitting Changes](#submitting-changes)
9. [Review Process](#review-process)
10. [Community](#community)

## Code of Conduct

The Phytovenomics project is committed to fostering an open and welcoming environment. Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

## Getting Started

Before you begin:

1. Make sure you have a GitHub account
2. Familiarize yourself with Git/GitHub workflow
3. Read the project documentation
4. Check the [issues](https://github.com/your-organization/phytovenomics/issues) for open tasks or create a new issue to discuss your proposed changes

## Development Environment

### Setting Up Your Environment

1. Clone the repository:
   ```bash
   git clone https://github.com/your-organization/phytovenomics.git
   cd phytovenomics
   ```

2. Set up a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install development dependencies:
   ```bash
   pip install -e ".[dev]"
   ```

4. Install pre-commit hooks:
   ```bash
   pre-commit install
   ```

### Installing Required Models

Some modules require pre-trained models that are not included in the repository due to their size. You can download them using our utility script:

```bash
python scripts/download_models.py --all
```

Or selectively download specific models:

```bash
python scripts/download_models.py --esm --igfold
```

## Branching Strategy

We follow a simplified Git flow approach:

- `main`: Production-ready code
- `develop`: Integration branch for new features
- `feature/*`: Individual feature branches
- `fix/*`: Bug fix branches
- `release/*`: Release preparation branches

Always create your working branches from `develop` and submit pull requests back to `develop`.

## Coding Standards

We follow standard Python conventions with some specific requirements:

### Style Guide

- Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) style guidelines
- Use [Black](https://black.readthedocs.io/) for code formatting
- Use [isort](https://pycqa.github.io/isort/) for import sorting
- Use [flake8](https://flake8.pycqa.org/) for linting

### Code Organization

- Keep modules focused on a single responsibility
- Use clear, descriptive names for functions, classes, and variables
- Document public APIs using docstrings (NumPy or Google style)
- Include type hints for function parameters and return values

### Example

```python
from typing import List, Optional

def calculate_binding_energy(antibody_sequence: str, toxin_sequence: str) -> float:
    """Calculate the binding energy between antibody and toxin sequences.
    
    Parameters
    ----------
    antibody_sequence : str
        The amino acid sequence of the antibody
    toxin_sequence : str
        The amino acid sequence of the toxin
        
    Returns
    -------
    float
        The estimated binding energy in kcal/mol
    """
    # Implementation...
    return energy_value
```

## Documentation

Documentation is crucial to our project. Please follow these guidelines:

1. **Code Documentation**:
   - Document all public modules, functions, classes, and methods
   - Update docstrings when changing function signatures
   - Include usage examples where appropriate

2. **README Files**:
   - Module-specific README files should explain the purpose, features, and usage of the module
   - Update READMEs when adding or changing functionality

3. **Tutorials and Examples**:
   - Consider adding tutorials for complex features
   - Update example code to reflect API changes

## Testing

All code contributions should include tests:

1. Create unit tests for new functions and methods
2. Ensure all tests pass before submitting changes
3. Aim for high test coverage, especially for critical functions

To run tests:

```bash
pytest
```

For coverage report:

```bash
pytest --cov=phytovenomics tests/
```

## Submitting Changes

1. Create a new branch from `develop`:
   ```bash
   git checkout develop
   git pull
   git checkout -b feature/your-feature-name
   ```

2. Make your changes, following our coding standards

3. Add tests for your changes

4. Run the full test suite to ensure nothing is broken:
   ```bash
   pytest
   ```

5. Commit your changes with descriptive commit messages:
   ```bash
   git commit -m "Add feature X to module Y"
   ```

6. Push your branch to GitHub:
   ```bash
   git push origin feature/your-feature-name
   ```

7. Create a pull request to the `develop` branch

## Review Process

All contributions go through review before being merged:

1. Automated checks must pass (tests, linting, etc.)
2. At least one maintainer must approve the changes
3. Address any feedback or requests for changes
4. Once approved, a maintainer will merge your contribution

## Community

Join our community channels for discussion and support:

- [GitHub Discussions](https://github.com/your-organization/phytovenomics/discussions)
- [Slack Channel](#) (Invitation available upon request)
- [Mailing List](#)

Thank you for contributing to Phytovenomics!