# Changelog
All notable changes to the Phytovenomics project will be documented in this file.

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - Development Version

### Added
- Proof-of-concept pipeline connecting all modules (venom intelligence, epitope discovery, cocktail formulator, antibody generator)
- Integration framework for modular toxin analysis
- Comprehensive setup system with environment validation
- Automated model downloading with integrity verification
- Configurable pipeline parameters via YAML configuration

### Changed
- Enhanced data preprocessing for ML model compatibility
- Improved epitope ranking algorithm for broader toxin coverage
- Refined antibody structure prediction workflow

### Fixed
- Memory leaks in large sequence processing
- Stability issues when handling toxin variant data

## [0.1.0] - Initial Release (Planned)

### Features
- Core venom intelligence module for toxin sequence analysis
- Epitope discovery pipeline for identifying immunogenic regions
- Basic cocktail formulator for antivenom design
- Antibody generator with structural validation
- Command-line interface for end-to-end pipeline execution
- Data preprocessing utilities for various input formats
- Documentation covering module design and usage
- Test suite for individual components and integration

### Dependencies
- Biopython 1.80+
- PyTorch 2.0+
- IgFold (custom package for immunoglobulin structure prediction)
- ESM-2 (embedded transformer models for toxin sequence analysis)
- NumPy, Pandas, and SciPy for data handling
- Matplotlib and Seaborn for visualization

[Unreleased]: https://github.com/yourusername/phytovenomics/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/yourusername/phytovenomics/releases/tag/v0.1.0