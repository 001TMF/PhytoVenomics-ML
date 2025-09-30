# Data Directory

This directory contains all input data, output results, and model weights for the antibody design pipeline.

## Directory Structure

```
data/
├── input/          # Input structures
│   ├── antigens/   # Snake toxin and target structures
│   ├── frameworks/ # Antibody framework templates
│   └── complexes/  # Pre-docked antibody-antigen complexes
│
├── output/         # Pipeline outputs (one directory per experiment)
│
└── models/         # Model weights (not tracked in git)
```

## Input Directory (`input/`)

### `antigens/`
Place your target structures here (snake toxins, other antigens).

**Supported formats**: PDB

**Example files**:
- `alpha_bungarotoxin.pdb`
- `taipoxin.pdb`
- `snake_toxin_7M6C.pdb`

### `frameworks/`
Place antibody framework templates here for design.

**Supported formats**: PDB (must contain variable regions)

**Example files**:
- `human_igg1_framework.pdb`
- `mouse_framework.pdb`

**Requirements**:
- Must contain heavy chain (H) and/or light chain (L)
- Should be numbered in any standard scheme (will be renumbered to Chothia)

### `complexes/`
Place pre-docked antibody-antigen complexes here.

**Supported formats**: PDB

**Example files**:
- `crystal_structure_7M6C.pdb`
- `validated_docked_pose.pdb`

**Requirements**:
- Must contain both antibody chains (H, L) and antigen chain(s)
- Chain IDs must be specified in config: `antibody_chains: [H, L]`, `antigen_chains: [A]`

## Output Directory (`output/`)

Pipeline creates one subdirectory per experiment:

```
output/
└── experiment_2024-01-15_cdr_h3/
    ├── prepared/           # Renumbered structures
    ├── docked/            # Docked complexes (Mode 1 only)
    ├── designs/           # Generated designs
    └── results.json       # Ranked candidates with metrics
```

**Note**: This directory is gitignored. Outputs are experiment-specific and can be large.

## Models Directory (`models/`)

Place downloaded model weights here.

**Example files**:
- `rfdiffusion_weights.pt`
- `iglm_weights.bin`

**Note**: Large model files are gitignored. Download separately as needed.

## Usage

### Mode 1: Separate Antigen + Framework

```bash
python -m antibody_pipeline \
  --antigen data/input/antigens/snake_toxin.pdb \
  --framework data/input/frameworks/human_framework.pdb \
  --output data/output/my_experiment
```

### Mode 2: Pre-Docked Complex

```bash
python -m antibody_pipeline \
  --complex data/input/complexes/crystal_structure.pdb \
  --antibody-chains H L \
  --antigen-chains A \
  --output data/output/my_experiment
```

### Mode 3: De Novo Design

```bash
python -m antibody_pipeline \
  --antigen data/input/antigens/snake_toxin.pdb \
  --output data/output/my_experiment
```

## Data Management

### Best Practices

1. **Organize by project**: Group related structures in subdirectories
2. **Use descriptive names**: `snake_toxin_species_date.pdb` instead of `file1.pdb`
3. **Track metadata**: Keep notes on structure sources and preparation
4. **Backup outputs**: Results can be large but may be valuable

### Example Organization

```
data/input/antigens/
├── snake_toxins/
│   ├── alpha_bungarotoxin/
│   │   ├── 1ABT.pdb
│   │   └── notes.txt
│   ├── taipoxin/
│   │   ├── 1TXN.pdb
│   │   └── notes.txt
│   └── ...
└── other_targets/
    └── ...
```

## Notes

- Input files should be clean PDB files
- Remove waters and ligands unless needed
- Check chain IDs before using pre-docked mode
- Output directories can grow large (GB+)