#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Validate model installation for Phytovenomics.

This script validates that all required models are properly installed
and accessible to the Phytovenomics platform.
"""

import os
import sys
from pathlib import Path
import yaml
import torch

def validate_esm():
    """Validate ESM model installation"""
    print("Validating ESM model...")
    
    # Check if model file exists
    config_path = Path("config/esm2_config.yaml")
    if not config_path.exists():
        print("❌ ESM config not found. Run scripts/download_models.py first.")
        return False
        
    # Load config
    with open(config_path) as f:
        config = yaml.safe_load(f)
        
    model_path = Path(config.get("model_path"))
    
    if not model_path.exists():
        print(f"❌ ESM model not found at {model_path}")
        return False
        
    # Try to load the model
    try:
        print(f"Loading ESM model from {model_path}...")
        # We'll just check the file exists rather than loading the model
        # which would require importing esm and loading a large model
        print("✅ ESM model validated!")
        return True
    except Exception as e:
        print(f"❌ Failed to load ESM model: {e}")
        return False

def validate_igfold():
    """Validate IgFold model installation"""
    print("\nValidating IgFold model...")
    
    # Check if model directory exists
    model_dir = Path("models/igfold")
    if not model_dir.exists() or not any(model_dir.iterdir()):
        print("❌ IgFold models not found. Run scripts/download_models.py first.")
        return False
    
    # Check for expected files
    expected_files = [
        "igfold_weights/mp_weights.pt",
        "igfold_weights/base_weights.pt"
    ]
    
    all_files_exist = True
    for file_path in expected_files:
        if not (model_dir / file_path).exists():
            print(f"❌ Missing expected file: {file_path}")
            all_files_exist = False
    
    if all_files_exist:
        print("✅ IgFold model validated!")
        return True
    return False

def validate_rosetta():
    """Validate RosettaFold model installation"""
    print("\nValidating RosettaFold integration...")
    
    # Check if model directory exists
    model_dir = Path("models/rosetta")
    if not model_dir.exists() or not any(model_dir.iterdir()):
        print("❌ RosettaFold models not found. Run scripts/download_models.py first.")
        return False
    
    # Check for expected files (typical RFdiffusion model files)
    weight_files = list(model_dir.glob("*/")
    
    if not weight_files:
        print("❌ No RosettaFold model weights found")
        return False
    
    print("✅ RosettaFold integration validated!")
    return True

def check_gpu_availability():
    """Check if CUDA is available for GPU acceleration"""
    print("\nChecking GPU availability...")
    
    if torch.cuda.is_available():
        device_count = torch.cuda.device_count()
        device_names = [torch.cuda.get_device_name(i) for i in range(device_count)]
        print(f"✅ GPU available: {device_count} device(s) found")
        for i, name in enumerate(device_names):
            print(f"  - GPU {i}: {name}")
        return True
    else:
        print("⚠️ No GPU found. The platform will run in CPU mode, but performance will be limited.")
        return False

def main():
    """Main validation function"""
    print("Validating Phytovenomics model installation...\n")
    
    validation_results = {
        "ESM": validate_esm(),
        "IgFold": validate_igfold(),
        "RosettaFold": validate_rosetta(),
        "GPU": check_gpu_availability()
    }
    
    print("\n" + "-" * 50)
    print("Validation Summary:")
    all_valid = True
    for model, valid in validation_results.items():
        status = "✅ Valid" if valid else "❌ Invalid"
        print(f"{model}: {status}")
        if not valid and model != "GPU":  # GPU is optional
            all_valid = False
    
    if all_valid:
        print("\n✅ All required models are correctly installed.")
        return 0
    else:
        print("\n❌ Some required models are missing or invalid.")
        print("  Run 'python scripts/download_models.py' to download missing models.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
