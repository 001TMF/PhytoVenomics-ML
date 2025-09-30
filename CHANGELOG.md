# Changelog

All notable changes to the Integrated Antibody Design Pipeline.

## [1.0.0] - 2024

### Added - Complete Integrated Pipeline

#### Core Features
- **Three input modes** for maximum flexibility:
  1. Separate antigen + framework (with docking)
  2. **Pre-docked antibody-antigen complex** (NEW - use existing crystal structures)
  3. Antigen only (de novo design)

#### Modules Implemented

**1. Preparation Module** (from DiffAb)
- Chothia antibody numbering using abnumber
- CDR identification and labeling (H1, H2, H3, L1, L2, L3)
- BioPython-based structure parsing
- Heavy atom extraction
- **Complex splitting for pre-docked structures** (NEW)

**2. Docking Module** (from DiffAb)
- HDock integration for protein-protein docking
- Epitope-guided docking
- Antibody-specific CDR H3 binding site
- **Skip docking option for pre-docked complexes** (NEW)

**3. Design Module** (RFdiffusion integration)
- RFdiffusion v3 wrapper
- Multiple design modes: CDR H3, single CDR, multiple CDRs, full
- Epitope conditioning
- Structure optimization

**4. Optimization Module** (from Germinal)
- IgLM-based sequence optimization
- Human-likeness scoring
- CDR-specific evaluation
- Mutation suggestions

**5. Filtering Module** (from Germinal)
- Structure prediction (Chai-1/AlphaFold3)
- Clash detection
- Interface metrics (contacts, BSA, pDockQ)
- Developability assessment (SAP score, PTM sites, charge patches)

**6. Pipeline Orchestrator**
- End-to-end workflow automation
- Configurable thresholds
- Multi-stage filtering
- Result ranking and export

**7. CLI and Configuration**
- Command-line interface with argparse
- YAML/JSON configuration support
- **Support for all three input modes** (NEW)
- Example configurations

#### Directory Structure
- Organized repository with clear hierarchy
- `antibody_pipeline/` - Main package
- `bin/` - External binaries
- `data/` - Input, output, and models
- `examples/` - Example configs and scripts
- `docs/` - Documentation
- `tests/` - Test suite (framework)
- `scripts/` - Utility scripts
- `logs/` - Log files

#### Documentation
- Complete README with all three modes
- Quick start guide
- Pipeline implementation details
- Integration notes
- Project structure documentation
- Example configurations for all modes

#### Installation
- `setup.py` for pip installation
- `requirements.txt` for dependencies
- `.gitignore` for version control
- Installation instructions for external tools

### Integration Strategy

**From DiffAb:**
- Structure preparation and renumbering
- CDR identification (Chothia scheme)
- HDock docking integration
- Structure parsing utilities

**From Germinal:**
- IgLM wrapper with gradient computation
- Filtering pipeline architecture
- Structure prediction integration
- Developability metrics

**New Implementation:**
- RFdiffusion v3 wrapper
- **Pre-docked complex handling**
- Complete pipeline orchestrator
- Unified configuration system

### Files Created

**Core Package** (~3,500 lines):
- `antibody_pipeline/pipeline.py` (450+ lines)
- `antibody_pipeline/cli.py` (280+ lines)
- `antibody_pipeline/preparation/` (4 files, 900+ lines)
- `antibody_pipeline/docking/` (2 files, 310+ lines)
- `antibody_pipeline/design/` (1 file, 330+ lines)
- `antibody_pipeline/optimization/` (1 file, 270+ lines)
- `antibody_pipeline/filtering/` (3 files, 820+ lines)
- `antibody_pipeline/config/` (1 file, 110+ lines)

**Documentation** (~2,500 lines):
- `README.md` - Main documentation with all modes
- `docs/QUICK_START.md` - 5-minute tutorial
- `docs/PIPELINE_SUMMARY.md` - Implementation details
- `docs/INTEGRATION_NOTES.md` - Integration strategy
- `PROJECT_STRUCTURE.md` - Directory organization
- `CHANGELOG.md` - This file

**Examples**:
- `examples/config_example.yaml` - Standard mode
- `examples/config_predocked.yaml` - **Pre-docked mode** (NEW)
- `examples/run_example.py` - Python API

**Configuration**:
- `setup.py` - Package installation
- `requirements.txt` - Dependencies
- `.gitignore` - Version control

### Usage Examples

#### Standard Mode (Antigen + Framework)
```bash
python -m antibody_pipeline \
  --antigen toxin.pdb \
  --framework template.pdb \
  --num-designs 10
```

#### Pre-Docked Mode (NEW)
```bash
python -m antibody_pipeline \
  --complex crystal.pdb \
  --antibody-chains H L \
  --antigen-chains A \
  --num-designs 10
```

#### De Novo Mode
```bash
python -m antibody_pipeline \
  --antigen toxin.pdb \
  --mode full \
  --num-designs 20
```

### Technical Details

**Languages**: Python 3.8+
**Key Dependencies**: PyTorch, BioPython, abnumber, iglm
**External Tools**: RFdiffusion, HDock (optional), Chai-1 (optional)

**Code Statistics**:
- 25+ Python files
- 3,500+ lines of production code
- 2,500+ lines of documentation
- 7 major modules

### Known Limitations

1. Requires separate RFdiffusion installation
2. HDock requires manual download
3. GPU recommended for most features
4. PyRosetta optional (requires license)
5. Chai/AF3 optional for structure prediction

### Future Enhancements

Potential additions for future versions:
- Additional structure predictors (ESMFold, OmegaFold)
- Support for single-domain antibodies (VHH)
- Batch processing for multiple targets
- Web interface
- Experimental validation integration
- Database storage for designs

## Version History

### [1.0.0] - 2024
- Initial release
- Complete integrated pipeline
- Three input modes including pre-docked complex support
- Full documentation and examples

---

For detailed information, see:
- [README.md](README.md) - Main documentation
- [docs/PIPELINE_SUMMARY.md](docs/PIPELINE_SUMMARY.md) - Implementation details
- [docs/INTEGRATION_NOTES.md](docs/INTEGRATION_NOTES.md) - Integration strategy