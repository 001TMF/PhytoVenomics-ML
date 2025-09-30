# Installation Guide

Complete installation instructions for the Integrated Antibody Design Pipeline.

## System Requirements

### Minimum Requirements
- **OS**: Linux, macOS, or Windows (WSL2)
- **Python**: ≥ 3.8
- **RAM**: 16 GB (32 GB recommended)
- **Storage**: 50 GB free space
- **GPU**: NVIDIA GPU with CUDA support (8+ GB VRAM recommended)

### Recommended Setup
- **OS**: Ubuntu 20.04+ or macOS 12+
- **Python**: 3.10
- **RAM**: 32+ GB
- **GPU**: NVIDIA A100, RTX 3090, or RTX 4090
- **CUDA**: 11.8 or 12.1

## Quick Installation (5 minutes)

```bash
# 1. Clone repository
git clone <repository-url>
cd Phytovenomics-ML

# 2. Create conda environment (recommended)
conda create -n antibody-design python=3.10
conda activate antibody-design

# 3. Install PyTorch (adjust CUDA version)
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia

# 4. Install pipeline
pip install -e .

# 5. Install RFdiffusion
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion
pip install -e .
cd ..

# 6. Validate installation
python scripts/validate_installation.py
```

## Detailed Installation

### Step 1: Python Environment

#### Option A: Conda (Recommended)

```bash
# Create environment
conda create -n antibody-design python=3.10
conda activate antibody-design

# Verify Python version
python --version  # Should be 3.10.x
```

#### Option B: venv

```bash
# Create virtual environment
python3.10 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Verify
python --version
```

### Step 2: PyTorch Installation

Choose the command for your system from https://pytorch.org/

#### CUDA 11.8
```bash
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
```

#### CUDA 12.1
```bash
conda install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia
```

#### CPU Only (Not Recommended)
```bash
conda install pytorch torchvision torchaudio cpuonly -c pytorch
```

#### Verify PyTorch
```python
python -c "import torch; print(f'PyTorch: {torch.__version__}'); print(f'CUDA available: {torch.cuda.is_available()}')"
```

### Step 3: Install Pipeline Package

```bash
cd Phytovenomics-ML

# Install in editable mode
pip install -e .

# Verify
python -c "from antibody_pipeline import AntibodyDesignPipeline; print('Success!')"
```

### Step 4: Install RFdiffusion

```bash
# Clone RFdiffusion
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion

# Install dependencies
pip install -e .

# Download weights (if needed)
# Follow instructions in RFdiffusion README

cd ..
```

### Step 5: Optional Dependencies

#### HDock (for Mode 1: separate antigen + framework)

```bash
# 1. Download from http://huanglab.phys.hust.edu.cn/software/hdocklite/
# 2. Extract archive
# 3. Copy binaries
cp /path/to/hdock bin/hdock
cp /path/to/createpl bin/createpl

# 4. Make executable
chmod +x bin/hdock
chmod +x bin/createpl

# 5. Test
./bin/hdock -h
```

#### Chai-1 (for structure prediction)

```bash
pip install chai_lab
```

#### PyRosetta (for relaxation and analysis)

PyRosetta requires a license. Get one from https://www.pyrosetta.org/downloads

```bash
# After obtaining license
pip install pyrosetta
```

### Step 6: Validate Installation

```bash
python scripts/validate_installation.py
```

Expected output:
```
============================================================
ANTIBODY DESIGN PIPELINE - INSTALLATION VALIDATION
============================================================
Checking Python version...
✓ Python 3.10.x

Checking Python packages...
✓ torch
✓ numpy
✓ biopython
✓ abnumber
✓ iglm
✓ transformers
✓ pyyaml

Checking antibody_pipeline installation...
✓ antibody_pipeline installed

Checking RFdiffusion...
✓ RFdiffusion found at ./RFdiffusion

Checking directory structure...
✓ antibody_pipeline/
✓ examples/
✓ docs/
...

============================================================
✓ Installation is complete and ready to use!
============================================================
```

## Troubleshooting

### Issue: CUDA not available

**Symptoms**: `torch.cuda.is_available()` returns `False`

**Solutions**:
1. Check NVIDIA driver: `nvidia-smi`
2. Reinstall PyTorch with correct CUDA version
3. Verify CUDA toolkit installation: `nvcc --version`

### Issue: Import errors

**Symptoms**: `ModuleNotFoundError` when importing

**Solutions**:
```bash
# Reinstall package
pip install -e . --force-reinstall

# Check installation
pip list | grep antibody
```

### Issue: HDock not found

**Symptoms**: "HDock binary not found"

**Solutions**:
1. Download HDock from official website
2. Place in `bin/` directory
3. Make executable: `chmod +x bin/hdock bin/createpl`
4. Or disable docking: `use_docking: false` in config

### Issue: Out of memory

**Symptoms**: CUDA out of memory errors

**Solutions**:
1. Reduce `num_designs` in config
2. Use CPU (slower): `device: cpu`
3. Use smaller batch sizes
4. Close other GPU applications

### Issue: RFdiffusion not found

**Symptoms**: "RFdiffusion not found at ./RFdiffusion"

**Solutions**:
```bash
# Clone RFdiffusion
git clone https://github.com/RosettaCommons/RFdiffusion.git

# Install
cd RFdiffusion && pip install -e . && cd ..

# Or specify path in config
rfdiffusion:
  path: /custom/path/to/RFdiffusion
```

## Platform-Specific Notes

### macOS

- HDock may require Rosetta 2 on Apple Silicon
- CUDA not available - use CPU mode or run on Linux
- Some features may be slower

### Windows (WSL2)

- Install WSL2 with Ubuntu 20.04+
- Follow Linux installation instructions
- CUDA should work if properly configured in WSL2

### Linux

- Most compatible platform
- CUDA support generally works out of the box
- Recommended for production use

## Verification Tests

### Test 1: Import Test
```python
from antibody_pipeline import AntibodyDesignPipeline, PipelineConfig
print("✓ Imports successful")
```

### Test 2: Configuration Test
```bash
python -m antibody_pipeline --create-config test_config.yaml
cat test_config.yaml
```

### Test 3: CLI Test
```bash
python -m antibody_pipeline --help
```

### Test 4: Module Tests
```python
from antibody_pipeline.preparation import AntibodyRenumberer
from antibody_pipeline.optimization import IgLMOptimizer
print("✓ All modules imported")
```

## Updating

### Update Pipeline
```bash
cd Phytovenomics-ML
git pull
pip install -e . --upgrade
```

### Update RFdiffusion
```bash
cd RFdiffusion
git pull
pip install -e . --upgrade
cd ..
```

### Update Dependencies
```bash
pip install -r requirements.txt --upgrade
```

## Uninstallation

```bash
# Remove package
pip uninstall antibody-design-pipeline

# Remove environment
conda env remove -n antibody-design

# Remove repository
rm -rf Phytovenomics-ML
```

## Getting Help

If installation fails:

1. Check this guide thoroughly
2. Run validation script: `python scripts/validate_installation.py`
3. Check [README.md](README.md) and [docs/QUICK_START.md](docs/QUICK_START.md)
4. Open GitHub issue with:
   - Error message
   - Output of validation script
   - OS and Python version
   - GPU details (if applicable)

## Next Steps

After successful installation:

1. **Read Quick Start**: `docs/QUICK_START.md`
2. **Try Examples**: `python examples/run_example.py`
3. **Create Config**: `python -m antibody_pipeline --create-config my_config.yaml`
4. **Run Pipeline**: `python -m antibody_pipeline --config my_config.yaml`

## Estimated Installation Time

- **Quick Install**: 5-10 minutes (without optional deps)
- **Full Install**: 15-30 minutes (with all optional deps)
- **Download Times**: Depends on internet speed
  - RFdiffusion: ~5 minutes
  - Model weights: 10-30 minutes (if needed)

## Disk Space Requirements

- **Pipeline code**: ~10 MB
- **RFdiffusion**: ~2 GB
- **Model weights**: ~5-10 GB
- **Example data**: ~100 MB
- **Working space**: 20+ GB (for outputs)
- **Total**: ~50 GB recommended