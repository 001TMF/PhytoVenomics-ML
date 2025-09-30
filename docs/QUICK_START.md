# Quick Start Guide

## Installation (5 minutes)

```bash
# 1. Install the pipeline
cd Phytovenomics-ML
pip install -e .

# 2. Install RFdiffusion
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion
pip install -e .
cd ..

# 3. Download HDock (optional, for docking)
# Get from: http://huanglab.phys.hust.edu.cn/software/hdocklite/
mkdir -p bin
# Place hdock and createpl binaries in bin/

# 4. Install Chai (optional, for structure prediction)
pip install chai_lab
```

## Minimal Example (2 minutes)

```python
from antibody_pipeline import AntibodyDesignPipeline, PipelineConfig

# Configure
config = PipelineConfig(
    antigen_pdb="your_toxin.pdb",
    framework_pdb="your_framework.pdb",  # Or None for de novo
    output_dir="./results",
    num_designs=5,
    design_mode="cdr_h3",
    rfdiffusion_path="./RFdiffusion",
)

# Run
pipeline = AntibodyDesignPipeline(config)
candidates = pipeline.run()

# Get best design
best = candidates[0]
print(f"Sequence: {best.sequence}")
print(f"PDB: {best.pdb_path}")
```

## Command Line (1 minute)

```bash
# Create config
python -m antibody_pipeline --create-config config.yaml

# Edit config.yaml with your paths

# Run
python -m antibody_pipeline --config config.yaml
```

## Configuration Template

```yaml
antigen_pdb: path/to/toxin.pdb
framework_pdb: path/to/framework.pdb
output_dir: ./results

design:
  mode: cdr_h3
  num_designs: 10

rfdiffusion:
  path: ./RFdiffusion

iglm:
  enabled: true
  species: '[HUMAN]'

device: cuda:0
```

## Output

Results saved to `output_dir/`:
```
results/
├── prepared/          # Renumbered structures
├── docked/           # Docked complexes
├── designs/          # Generated designs
└── results.json      # Ranked candidates
```

## Common Use Cases

### 1. Design CDR H3 for known epitope
```python
config = PipelineConfig(
    antigen_pdb="toxin.pdb",
    framework_pdb="framework.pdb",
    epitope_residues=["A45", "A48", "A52"],
    design_mode="cdr_h3",
    num_designs=20,
)
```

### 2. De novo design (no framework)
```python
config = PipelineConfig(
    antigen_pdb="toxin.pdb",
    framework_pdb=None,  # Design from scratch
    design_mode="full",
    num_designs=50,
)
```

### 3. Fast screening (no IgLM)
```python
config = PipelineConfig(
    antigen_pdb="toxin.pdb",
    framework_pdb="framework.pdb",
    use_iglm_optimization=False,  # Skip IgLM
    num_designs=100,
)
```

## Troubleshooting

### "RFdiffusion not found"
- Ensure RFdiffusion is installed and `rfdiffusion_path` is correct
- Check: `ls RFdiffusion/scripts/run_inference.py`

### "HDock binary not found"
- Download HDock from http://huanglab.phys.hust.edu.cn/software/hdocklite/
- Place in `./bin/` or specify path in config

### "CUDA out of memory"
- Reduce `num_designs`
- Use smaller batch sizes
- Or use CPU: `device: cpu` (slower)

### "No designs passed filters"
- Check thresholds in config (may be too strict)
- Review design quality in `designs/` directory
- Adjust: `max_clashes`, `min_plddt`, `max_sap_score`

## Next Steps

1. **Read full docs**: See [README_PIPELINE.md](README_PIPELINE.md)
2. **Check examples**: See [examples/](examples/)
3. **View architecture**: See [PIPELINE_SUMMARY.md](PIPELINE_SUMMARY.md)
4. **Integration details**: See [INTEGRATION_NOTES.md](INTEGRATION_NOTES.md)

## Getting Help

- Check documentation files in repo
- Review example scripts in `examples/`
- Open GitHub issue for bugs