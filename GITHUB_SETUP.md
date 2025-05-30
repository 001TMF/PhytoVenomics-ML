# Phytovenomics GitHub Setup Guide

## Repository Setup

This guide will walk you through setting up and pushing the Phytovenomics codebase to GitHub.

### 1. Create a New GitHub Repository

1. Go to [GitHub](https://github.com/) and log in to your account
2. Click the "+" icon in the top right corner and select "New repository"
3. Name your repository: `phytovenomics`
4. Add a brief description: "Machine learning platform for plant-based antivenom development"
5. Choose visibility (Public or Private)
6. **DO NOT** initialize the repository with README, .gitignore, or license files (we already have these)
7. Click "Create repository"

### 2. Prepare the Local Repository

The repository has already been initialized with Git. Before pushing, ensure that all generated files and directories that shouldn't be tracked are properly excluded:

```bash
# Review files that would be committed
git status

# Ensure large data files are properly excluded
find . -type f -size +10M | grep -v ".git"
```

Our `.gitignore` file should handle most exclusions. If you notice any additional files that should be excluded, add them to the `.gitignore` file.

### 3. Configure Git

Set your Git user information if you haven't already:

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

### 4. Stage and Commit Files

```bash
# Add all files to staging
git add .

# Commit the changes
git commit -m "Initial commit of Phytovenomics platform"
```

### 5. Connect to the Remote Repository

Replace `YOUR_USERNAME` with your GitHub username:

```bash
git remote add origin https://github.com/YOUR_USERNAME/phytovenomics.git
git branch -M main
```

### 6. Push the Repository

```bash
git push -u origin main
```

You may need to authenticate with your GitHub credentials or a Personal Access Token.

## Repository Structure

The Phytovenomics repository has the following key components:

### Core Components

- **phytovenomics/**: Main package directory
  - **bin/**: Command-line interface entry points
  - **setup/**: Setup and configuration utilities
  - **cli/**: Command-line interface implementation

### Documentation

- **docs/**: Comprehensive documentation
  - **modules/**: Documentation for individual modules
  - **technical_integration_report.md**: Technical report on module integration

### Data

- **data/**: Data files used throughout the project
  - **snake_toxins/**: Toxin sequence and structure data
  - **antibody_structures/**: Antibody structural data
  - **toxin_antibody_binding/**: Binding affinity data
  - **ml_ready/**: Processed datasets ready for ML training
  - **reports/**: Analysis reports including ML data alignment
  - **visualizations/**: Generated visualizations and plots

### Development

- **tests/**: Test files including the proof-of-concept pipeline test
- **.github/workflows/**: CI/CD configuration files
- **setup.py**: Package installation script
- **requirements.txt**: Package dependencies

### Source Code References

- **src/**: Additional source code modules
  - **structure_prediction/**: Hybrid structure prediction implementation
- **python_template/**: Template code and demonstrations

## Key Files

- **README.md**: Project overview and introduction
- **README_POC.md**: Proof-of-concept pipeline documentation
- **phytovenomics_system_design.md**: System architecture documentation
- **phytovenomics_class_diagram.mermaid**: Class structure diagram
- **phytovenomics_sequence_diagram.mermaid**: Process flow diagram

## Getting Started After Clone

After cloning the repository, set up the project by running:

```bash
# Install the package in development mode
pip install -e .[dev]

# Run the proof-of-concept test
python tests/proof_of_concept_test.py
```

## CI/CD

The repository includes GitHub Actions workflows for continuous integration. These workflows automatically run tests and linting checks when changes are pushed to the repository.
