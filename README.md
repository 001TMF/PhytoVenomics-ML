# Integrated Antibody Design Pipeline

Complete end-to-end antibody design system for snake toxin neutralization, combining DiffAb preparation, HDock docking, RFdiffusion design, IgLM optimization, and comprehensive filtering.

## Features

### Three Input Modes

1. **Separate antigen + framework** (will dock using HDock)
2. **Pre-docked antibody-antigen complex** (use existing crystal structure or docked pose)
3. **Antigen only** (de novo antibody design)

### Core Capabilities

- **Antibody Preparation** (DiffAb): Chothia numbering, CDR identification, structure parsing
- **Flexible Docking** (HDock): New docking or use pre-docked complexes, epitope-guided
- **Structure-Based Design** (RFdiffusion): CDR H3, multiple CDRs, or full variable region design
- **Sequence Optimization** (IgLM/Germinal): Human-likeness scoring, CDR evaluation
- **Comprehensive Filtering** (Germinal): Structure prediction, interface metrics, developability

## Quick Start

### Installation

```bash
# 1. Install package
cd Phytovenomics-ML
pip install -e .

# 2. Install RFdiffusion
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion && pip install -e . && cd ..

# 3. Optional: HDock for docking
# Download from http://huanglab.phys.hust.edu.cn/software/hdocklite/
# Place binaries in bin/

# 4. Optional: Chai for structure prediction
pip install chai_lab
```

### Mode 1: Separate Antigen + Framework

```bash
python -m antibody_pipeline \
  --antigen data/input/snake_toxin.pdb \
  --framework data/input/human_framework.pdb \
  --mode cdr_h3 \
  --num-designs 10 \
  --output ./results
```

### Mode 2: Pre-Docked Complex ⭐ NEW

```bash
python -m antibody_pipeline \
  --complex data/input/crystal_structure.pdb \
  --antibody-chains H L \
  --antigen-chains A \
  --mode cdr_h3 \
  --num-designs 10 \
  --output ./results
```

Or use config file:
```bash
python -m antibody_pipeline --config examples/config_predocked.yaml
```

### Mode 3: De Novo Design

```bash
python -m antibody_pipeline \
  --antigen data/input/snake_toxin.pdb \
  --mode full \
  --num-designs 20 \
  --output ./results
```

## Python API

```python
from antibody_pipeline import AntibodyDesignPipeline, PipelineConfig

# Pre-docked complex example
config = PipelineConfig(
    complex_pdb="data/input/antibody_toxin_complex.pdb",
    antibody_chains=["H", "L"],
    antigen_chains=["A"],
    output_dir="./results",
    num_designs=10,
    design_mode="cdr_h3",
    use_docking=False,
    rfdiffusion_path="./RFdiffusion",
)

pipeline = AntibodyDesignPipeline(config)
candidates = pipeline.run()

# View best design
best = candidates[0]
print(f"Best: {best.sequence}")
print(f"pLDDT: {best.plddt}")
print(f"SAP: {best.sap_score}")
```

## Directory Structure

```
antibody_pipeline/          # Main package
├── preparation/            # Structure prep (DiffAb)
├── docking/                # HDock integration
├── design/                 # RFdiffusion wrapper
├── optimization/           # IgLM optimizer (Germinal)
├── filtering/              # Quality filters (Germinal)
└── pipeline.py             # Main orchestrator

bin/                        # External binaries (hdock, createpl)
data/                       # Data (input/, output/, models/)
examples/                   # Example configs and scripts
docs/                       # Documentation
```

## Output

Results in `output_dir/results.json`:

```json
{
  "candidates": [
    {
      "design_id": "design_003",
      "rank": 1,
      "sequence": "EVQLVES...",
      "pdb_path": "./results/designs/design_003.pdb",
      "cdr_sequences": {"H1": "...", "H2": "...", "H3": "..."},
      "iglm_log_likelihood": -45.2,
      "plddt": 82.5,
      "sap_score": 23.4,
      "passed_filters": true
    }
  ]
}
```

## Configuration

Create default config:
```bash
python -m antibody_pipeline --create-config my_config.yaml
```

Example configurations:
- `examples/config_example.yaml` - Standard mode (separate antigen + framework)
- `examples/config_predocked.yaml` - Pre-docked complex mode
- `examples/run_example.py` - Python API example

## Documentation

- **[Quick Start](docs/QUICK_START.md)** - 5-minute tutorial
- **[Pipeline Summary](docs/PIPELINE_SUMMARY.md)** - Implementation details
- **[Integration Notes](docs/INTEGRATION_NOTES.md)** - How DiffAb/Germinal/RFdiffusion are integrated

## Common Workflows

### Redesign CDR H3 of Crystal Structure

```bash
python -m antibody_pipeline \
  --complex crystal.pdb \
  --antibody-chains H L \
  --antigen-chains A \
  --mode cdr_h3 \
  --num-designs 20
```

### Design for Known Epitope

```bash
python -m antibody_pipeline \
  --antigen toxin.pdb \
  --framework template.pdb \
  --epitope A45 A48 A52 \
  --num-designs 50
```

### High-Throughput Screening

```bash
python -m antibody_pipeline \
  --antigen toxin.pdb \
  --framework framework.pdb \
  --num-designs 100 \
  --no-iglm  # Skip optimization for speed
```

## Requirements

- Python ≥ 3.8
- PyTorch ≥ 2.0
- BioPython ≥ 1.79
- abnumber ≥ 0.3.0
- iglm ≥ 0.4.0
- Optional: chai_lab, PyRosetta

## Citation

If you use this pipeline, please cite:
- **DiffAb**: Luo et al.
- **Germinal**: [Paper reference]
- **RFdiffusion**: Watson et al., Nature 2023
- **IgLM**: Shuai et al., Nature Comm 2023

## License

MIT License

## Support

- [Documentation](docs/)
- [Examples](examples/)
- [GitHub Issues](https://github.com/your-repo/issues)