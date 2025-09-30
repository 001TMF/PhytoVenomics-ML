# Project Structure

## Directory Organization

```
Phytovenomics-ML/
│
├── antibody_pipeline/              # Main Python package
│   ├── __init__.py
│   ├── __main__.py                 # Entry point (python -m antibody_pipeline)
│   ├── pipeline.py                 # Main orchestrator (~450 lines)
│   ├── cli.py                      # Command-line interface (~280 lines)
│   │
│   ├── preparation/                # Structure preparation (DiffAb-based)
│   │   ├── __init__.py
│   │   ├── constants.py            # CDR definitions, AA types (~180 lines)
│   │   ├── parser.py               # BioPython structure parsing (~270 lines)
│   │   ├── renumber.py             # Chothia renumbering (~230 lines)
│   │   └── cdr_identifier.py       # CDR identification (~230 lines)
│   │
│   ├── docking/                    # Docking module (DiffAb-based)
│   │   ├── __init__.py
│   │   ├── base.py                 # Abstract interface (~60 lines)
│   │   └── hdock.py                # HDock wrapper (~250 lines)
│   │
│   ├── design/                     # Design module (RFdiffusion)
│   │   ├── __init__.py
│   │   └── designer.py             # RFdiffusion wrapper (~330 lines)
│   │
│   ├── optimization/               # Optimization (Germinal-based)
│   │   ├── __init__.py
│   │   └── iglm_optimizer.py       # IgLM wrapper (~270 lines)
│   │
│   ├── filtering/                  # Filtering (Germinal-based)
│   │   ├── __init__.py
│   │   ├── structure_filters.py    # Structure prediction (~240 lines)
│   │   ├── interface_metrics.py    # Interface analysis (~300 lines)
│   │   └── developability.py       # Developability (~280 lines)
│   │
│   ├── config/                     # Configuration management
│   │   ├── __init__.py
│   │   └── parser.py               # YAML/JSON parser (~110 lines)
│   │
│   └── utils/                      # Utilities (if needed)
│       └── __init__.py
│
├── bin/                            # External binaries
│   ├── hdock                       # HDock executable (download separately)
│   └── createpl                    # HDock post-processing
│
├── data/                           # Data directory
│   ├── input/                      # Input structures
│   │   ├── antigens/               # Snake toxin structures
│   │   ├── frameworks/             # Antibody frameworks
│   │   └── complexes/              # Pre-docked complexes
│   ├── output/                     # Pipeline outputs
│   │   └── [experiment_name]/      # Per-experiment results
│   └── models/                     # Model weights
│       └── (large files not in git)
│
├── examples/                       # Example configurations
│   ├── config_example.yaml         # Standard mode (antigen + framework)
│   ├── config_predocked.yaml       # Pre-docked complex mode
│   └── run_example.py              # Python API example
│
├── docs/                           # Documentation
│   ├── QUICK_START.md              # 5-minute tutorial
│   ├── PIPELINE_SUMMARY.md         # Implementation details
│   └── INTEGRATION_NOTES.md        # Integration strategy
│
├── tests/                          # Tests (to be implemented)
│   ├── __init__.py
│   ├── test_preparation.py
│   ├── test_docking.py
│   ├── test_design.py
│   ├── test_optimization.py
│   ├── test_filtering.py
│   └── test_pipeline.py
│
├── scripts/                        # Utility scripts
│   ├── download_models.py          # Download model weights
│   └── validate_installation.py   # Check installation
│
├── logs/                           # Log files
│   └── *.log
│
├── .gitignore                      # Git ignore rules
├── setup.py                        # Package installation
├── requirements.txt                # Python dependencies
├── README.md                       # Main documentation
└── PROJECT_STRUCTURE.md            # This file
```

## External Dependencies (Not in Repository)

These are managed separately and referenced via configuration:

```
RFdiffusion/                        # Clone from GitHub
diffab-main/                        # Reference implementation (for comparison)
germinal-main/                      # Reference implementation (for comparison)
```

## Input/Output Flow

### Mode 1: Separate Antigen + Framework

```
data/input/antigens/toxin.pdb
data/input/frameworks/template.pdb
            ↓
    [Pipeline Processing]
            ↓
data/output/experiment1/
    ├── prepared/
    │   ├── framework_chothia.pdb
    │   └── ...
    ├── docked/
    │   ├── docked_0.pdb
    │   └── ...
    ├── designs/
    │   ├── design_000.pdb
    │   └── ...
    └── results.json
```

### Mode 2: Pre-Docked Complex

```
data/input/complexes/crystal.pdb
            ↓
    [Pipeline Processing]
            ↓
data/output/experiment2/
    ├── prepared/
    │   ├── antibody_from_complex.pdb
    │   ├── antigen_from_complex.pdb
    │   └── antibody_chothia.pdb
    ├── designs/
    │   ├── design_000.pdb
    │   └── ...
    └── results.json
```

## Code Organization

### Module Hierarchy

```
antibody_pipeline.pipeline.AntibodyDesignPipeline
    ├── uses: preparation.*
    ├── uses: docking.*
    ├── uses: design.*
    ├── uses: optimization.*
    └── uses: filtering.*

antibody_pipeline.cli
    └── calls: AntibodyDesignPipeline
```

### Key Classes

- `PipelineConfig` - Configuration dataclass
- `AntibodyDesignPipeline` - Main orchestrator
- `DesignCandidate` - Result container
- `AntibodyRenumberer` - Chothia renumbering
- `HDockAntibody` - Docking engine
- `AntibodyDesigner` - RFdiffusion wrapper
- `IgLMOptimizer` - Sequence optimizer
- `StructureFilter` - Structure prediction
- `InterfaceMetrics` - Interface analysis
- `DevelopabilityFilter` - Developability assessment

## File Naming Conventions

- **Python modules**: `lowercase_with_underscores.py`
- **Classes**: `PascalCase`
- **Functions**: `snake_case`
- **Constants**: `UPPER_CASE`
- **Config files**: `config_description.yaml`
- **PDB files**: `descriptive_name.pdb`
- **Output dirs**: `experiment_name_date/`

## Code Statistics

- **Total Python files**: 25+
- **Total lines of code**: ~3,500+
- **Main modules**: 7 (preparation, docking, design, optimization, filtering, config, pipeline)
- **Documentation**: ~2,000+ lines
- **Examples**: 3 files

## Development Workflow

```bash
# 1. Clone repository
git clone <repo-url>
cd Phytovenomics-ML

# 2. Install in development mode
pip install -e .

# 3. Install external dependencies
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion && pip install -e . && cd ..

# 4. Download HDock (if needed)
# Place binaries in bin/

# 5. Run tests (when implemented)
pytest tests/

# 6. Run example
python examples/run_example.py

# 7. Or use CLI
python -m antibody_pipeline --config examples/config_example.yaml
```

## Git Strategy

### Tracked Files

- Python source code (`antibody_pipeline/`)
- Documentation (`docs/`, `README.md`)
- Examples (`examples/`)
- Tests (`tests/`)
- Configuration (`setup.py`, `requirements.txt`, `.gitignore`)

### Ignored Files (see .gitignore)

- `__pycache__/`, `*.pyc`
- `.venv/`, `venv/`
- `data/output/*/` (outputs)
- `data/models/*.pt` (large weights)
- `logs/*.log`
- `.idea/`, `.vscode/`
- `diffab-main/`, `germinal-main/`, `RFdiffusion/` (external)
- `bin/hdock`, `bin/createpl` (binaries)

## Data Management

### Input Data Organization

```
data/input/
├── antigens/
│   ├── snake_toxins/
│   │   ├── alpha_bungarotoxin.pdb
│   │   ├── taipoxin.pdb
│   │   └── ...
│   └── other_targets/
├── frameworks/
│   ├── human_igg1_vh.pdb
│   ├── human_igg1_vl.pdb
│   └── ...
└── complexes/
    ├── 7M6C.pdb
    ├── crystal_structures/
    └── docked_poses/
```

### Output Data Organization

```
data/output/
├── experiment_2024-01-15_cdr_h3/
│   ├── prepared/
│   ├── docked/
│   ├── designs/
│   └── results.json
├── experiment_2024-01-16_full_design/
└── ...
```

## Maintenance

### Adding New Features

1. Create module in appropriate directory
2. Add tests in `tests/`
3. Update documentation
4. Add example to `examples/`
5. Update `README.md`

### Updating External Dependencies

1. Check DiffAb/Germinal/RFdiffusion releases
2. Update integration code if needed
3. Update `requirements.txt`
4. Run tests

### Version Control

- Use semantic versioning (MAJOR.MINOR.PATCH)
- Tag releases in git
- Update `setup.py` version
- Update `CHANGELOG.md` (when created)