#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Download required models for the Phytovenomics platform.

This script downloads and sets up all necessary models required by the Phytovenomics
platform, including ESMFold, IgFold, and RosettaFold integration components.

Usage:
    python download_models.py [--all] [--esm] [--igfold] [--rosetta]

Options:
    --all       Download all models (default)
    --esm       Download only ESM models
    --igfold    Download only IgFold models
    --rosetta   Download only RosettaFold models
"""

import argparse
import os
import sys
import urllib.request
import tarfile
import zipfile
from pathlib import Path
import subprocess
import shutil
from tqdm import tqdm

# Define model URLs and directories
MODEL_DIR = Path("models")
MODELS = {
    "esm": {
        "url": "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t36_3B_UR50D.pt",
        "filename": "esm2_t36_3B_UR50D.pt",
        "dir": MODEL_DIR / "esm"
    },
    "igfold": {
        "url": "https://github.com/Graylab/IgFold/releases/download/v1.0.0/igfold_weights.zip",
        "filename": "igfold_weights.zip",
        "dir": MODEL_DIR / "igfold"
    },
    "rosetta": {
        "url": "https://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/weights.tar.gz",
        "filename": "weights.tar.gz",
        "dir": MODEL_DIR / "rosetta"
    }
}

class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

def download_file(url, filename):
    """Download a file with a progress bar"""
    with DownloadProgressBar(unit='B', unit_scale=True, miniters=1, desc=filename) as t:
        urllib.request.urlretrieve(url, filename=filename, reporthook=t.update_to)

def extract_archive(filename, extract_dir):
    """Extract tar.gz or zip archive"""
    print(f"Extracting {filename} to {extract_dir}...")
    
    if filename.endswith(".tar.gz") or filename.endswith(".tgz"):
        with tarfile.open(filename, "r:gz") as tar:
            tar.extractall(path=extract_dir)
    elif filename.endswith(".zip"):
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
    else:
        print(f"Unsupported archive format: {filename}")
        return False
        
    return True

def setup_model(model_key):
    """Download and setup a specific model"""
    model_info = MODELS[model_key]
    model_dir = model_info["dir"]
    model_file = model_dir / model_info["filename"]
    
    # Create model directory if it doesn't exist
    model_dir.mkdir(parents=True, exist_ok=True)
    
    # Download model file if it doesn't exist
    if not model_file.exists():
        print(f"Downloading {model_key} model...")
        download_file(model_info["url"], model_file)
    else:
        print(f"{model_key} model already downloaded.")
    
    # Extract archives if needed
    if model_file.name.endswith((".tar.gz", ".tgz", ".zip")):
        extract_archive(model_file, model_dir)

def setup_environment():
    """Create necessary config files and environment setup"""
    # Create config directory if it doesn't exist
    config_dir = Path("config")
    config_dir.mkdir(exist_ok=True)
    
    # Update model paths in config files
    esm_config = config_dir / "esm2_config.yaml"
    if not esm_config.exists():
        with open(esm_config, "w") as f:
            f.write(f"model_path: {MODELS['esm']['dir'] / MODELS['esm']['filename']}\n")
            f.write("use_gpu: true\n")
    
    # Create model registry file
    with open(MODEL_DIR / "registry.txt", "w") as f:
        for model_key, model_info in MODELS.items():
            f.write(f"{model_key}: {model_info['dir']}\n")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Download required models for Phytovenomics")
    parser.add_argument("--all", action="store_true", help="Download all models (default)")
    parser.add_argument("--esm", action="store_true", help="Download only ESM models")
    parser.add_argument("--igfold", action="store_true", help="Download only IgFold models")
    parser.add_argument("--rosetta", action="store_true", help="Download only RosettaFold models")
    
    args = parser.parse_args()
    
    # If no specific model is selected, download all
    if not (args.esm or args.igfold or args.rosetta):
        args.all = True
    
    # Create the models directory if it doesn't exist
    MODEL_DIR.mkdir(exist_ok=True)
    
    # Download selected models
    if args.all or args.esm:
        setup_model("esm")
    
    if args.all or args.igfold:
        setup_model("igfold")
    
    if args.all or args.rosetta:
        setup_model("rosetta")
    
    # Setup environment and configuration
    setup_environment()
    
    print("\nModel download complete!")
    print("Run 'python -m phytovenomics.validate_models' to verify installation.")

if __name__ == "__main__":
    main()
