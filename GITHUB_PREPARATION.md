# Phytovenomics GitHub Preparation Guide

This document provides specific instructions for preparing the Phytovenomics codebase for GitHub, with particular attention to data files and version control best practices.

## Data Files Management

### Current .gitignore Configuration

The project's `.gitignore` file already excludes many large data files, including:
- CSV files (*.csv)
- Pickle files (*.pkl)
- HDF5 files (*.h5, *.hdf5)
- Log files (*.log)
- Temporary and build directories (build/, dist/, tmp/)

### Git LFS Recommendation

For the Phytovenomics project, we recommend using Git Large File Storage (Git LFS) for:

1. **Essential data files that need version control**:
   - Reference datasets in `data/ml_ready/` that are crucial for reproducibility
   - Pre-trained model weights (if under 100MB)
   - Benchmark datasets for system validation

2. **Files to exclude completely**:
   - Raw data files in `data/processed/` and `data/snake_toxins/` 
   - Large PDB structure files
   - Generated visualization outputs
   - Local test results

## Git LFS Setup

To set up Git LFS for the Phytovenomics project:

```bash
# Install Git LFS
# For Ubuntu/Debian
sudo apt-get install git-lfs

# For macOS (using Homebrew)
brew install git-lfs

# Initialize Git LFS
git lfs install

# Track specific file types with Git LFS
git lfs track "data/ml_ready/*.csv"
git lfs track "data/reports/*.png"
git lfs track "models/*.h5"

# Make sure .gitattributes is tracked
git add .gitattributes
```

## Repository Size Considerations

GitHub has the following limitations:
- Recommended repository size: < 1GB
- Maximum file size: 100MB (without Git LFS)
- Maximum Git LFS file size: 2GB

The Phytovenomics project contains:
- Code (~10MB)
- Documentation (~5MB)
- Data files (potentially >100MB)
- Model files (potentially >100MB)

## Preparing for GitHub

Run the provided `prepare_repo.sh` script to:
1. Clean Python cache files
2. Check for large files
3. Identify potential sensitive data
4. Ensure scripts are executable

```bash
# Run the preparation script
./prepare_repo.sh

# Initialize Git repository if needed
git init

# Stage all files
git add .

# Initial commit
git commit -m "Initial commit of Phytovenomics platform"
```


## Data Directories Strategy

For GitHub deployment, we recommend:

1. **Include in Repository (with Git LFS)**:
   - `data/ml_ready/*.json` (statistics and metadata files)
   - `data/reports/*.md` (report files)
   - `data/visualizations/*.png` (important visualizations under 5MB)

2. **Provide Download Instructions**:
   - Large datasets
   - PDB structure files
   - Full visualization sets

## Sensitive Information Check

Before pushing to GitHub, verify no sensitive information exists in:
- Configuration files
- Notebooks
- Test files
- Documentation

## GitHub Actions

The repository includes a GitHub Actions workflow in `.github/workflows/ci.yml` that will:
- Run tests on each push/pull request
- Check code formatting with black
- Verify package builds correctly

## Getting Started After Clone

Update the README.md to include:

1. Installation instructions:
   ```bash
   # Clone the repository
   git clone https://github.com/yourusername/phytovenomics.git
   cd phytovenomics
   
   # Install the package in development mode
   pip install -e .[dev]
   
   # Download required models (if applicable)
   python scripts/download_models.py
   ```

2. Instructions for running the proof-of-concept pipeline:
   ```bash
   # Run the proof-of-concept pipeline test
   python tests/proof_of_concept_test.py
   ```

3. Data download instructions for any large files excluded from the repository

## Final Checklist

- [ ] Cleaned repository using prepare_repo.sh
- [ ] Configured Git LFS for necessary large files
- [ ] Verified .gitignore excludes appropriate files
- [ ] Checked for and removed sensitive information
- [ ] Updated README.md with complete setup instructions
- [ ] Included LICENSE file
- [ ] Created GitHub repository
- [ ] Pushed code to GitHub
- [ ] Verified GitHub Actions CI workflow runs successfully