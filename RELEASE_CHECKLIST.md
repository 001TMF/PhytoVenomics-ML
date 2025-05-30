# Phytovenomics Release Checklist

This document outlines the process for creating new releases of the Phytovenomics project on GitHub.

## Version Numbering

Phytovenomics follows [Semantic Versioning](https://semver.org/) (SemVer):

- **MAJOR version** (X.0.0): Incompatible API changes
- **MINOR version** (0.X.0): Add functionality in a backward-compatible manner
- **PATCH version** (0.0.X): Backward-compatible bug fixes

## Pre-Release Checklist

Before creating a new release:

- [ ] Update version number in `phytovenomics/_version.py`
- [ ] Update CHANGELOG.md with all significant changes
- [ ] Run the full test suite: `pytest`
- [ ] Check code coverage: `pytest --cov=phytovenomics`
- [ ] Run linting checks: `flake8` and `black --check phytovenomics tests`
- [ ] Update documentation to reflect new features or changes
- [ ] Verify that the README.md is up-to-date
- [ ] Verify all CI/CD workflows pass on the main branch

## Creating a GitHub Release

1. **Tag the release**:
   ```bash
   git checkout main
   git pull
   # Create an annotated tag
   git tag -a v0.1.0 -m "Version 0.1.0: Initial release"
   # Push the tag
   git push origin v0.1.0
   ```

2. **Create a GitHub Release**:
   - Go to the repository on GitHub
   - Navigate to "Releases"
   - Click "Draft a new release"
   - Select the tag you just created
   - Add a title (e.g., "Phytovenomics v0.2.0")
   - Include detailed release notes from CHANGELOG.md
   - Attach any relevant binaries or assets

3. **Release notes should include**:
   - Summary of major features and improvements
   - Detailed list of changes, organized by type:
     - ✨ New Features
     - 🐛 Bug Fixes
     - 📚 Documentation
     - 🛠 Internal Improvements
   - Breaking changes (if any)
   - Migration guide for users (if needed)
   - Acknowledgments for contributors

## Post-Release Tasks

- [ ] Announce the release to relevant channels
- [ ] Update any demonstration or example code to use the latest version
- [ ] Create a new milestone for the next version
- [ ] Verify the package can be installed from GitHub
- [ ] Update documentation website (if applicable)

## Release Branches

For major releases (X.0.0), consider creating a release branch:

```bash
# Create a release branch from main
git checkout main
git pull
git checkout -b release/1.0.x

# Push the branch
git push -u origin release/1.0.x
```


## Hotfixes

For urgent fixes to a released version:

1. Create a hotfix branch from the release tag:
   ```bash
   git checkout v1.0.0
   git checkout -b hotfix/v1.0.1
   # Make necessary changes
   # Update version number to v1.0.1
   ```

2. Make necessary fixes, update version number

3. Submit a PR to both:
   - The release branch (if applicable)
   - The main branch

4. After approval, tag and release the hotfix

## Testing Before Release

Test the installable package before final release:

```bash
# Clean build directories
rm -rf build/ dist/ *.egg-info/

# Create a source distribution and wheel
python setup.py sdist bdist_wheel

# Create a clean test environment
python -m venv test_env
source test_env/bin/activate  # On Windows: test_env\Scripts\activate

# Install the package from the built wheel
pip install dist/phytovenomics-*.whl

# Run a basic import test
python -c "import phytovenomics; print(phytovenomics.__version__)"

# Run the proof-of-concept test
python -m phytovenomics.tests.proof_of_concept_test

# Deactivate the environment
deactivate
```


## Special Considerations for Phytovenomics

### Model Dependencies

- [ ] Ensure model download scripts point to the correct version of models
- [ ] Update model compatibility notes if needed
- [ ] Test the model download functionality on a clean system

### Data Files

- [ ] Review any changes to data formats
- [ ] Update documentation for data inputs/outputs
- [ ] Make sure example data files are compatible with the release

### Pipeline Integration

- [ ] Verify that the proof-of-concept pipeline runs end-to-end 
- [ ] Test with both minimal and full data configurations
- [ ] Check for any regressions in pipeline outputs