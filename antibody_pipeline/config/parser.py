"""
Configuration file parser

Supports YAML and JSON configuration files.
"""

import os
import json
import yaml
from typing import Dict, Any


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from YAML or JSON file.

    Args:
        config_path: Path to configuration file

    Returns:
        Configuration dictionary

    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If file format not supported
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    ext = os.path.splitext(config_path)[1].lower()

    if ext in ['.yaml', '.yml']:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    elif ext == '.json':
        with open(config_path, 'r') as f:
            config = json.load(f)
    else:
        raise ValueError(f"Unsupported config format: {ext}. Use .yaml, .yml, or .json")

    return config


def save_config(config: Dict[str, Any], output_path: str):
    """
    Save configuration to YAML or JSON file.

    Args:
        config: Configuration dictionary
        output_path: Path to save configuration

    Raises:
        ValueError: If file format not supported
    """
    ext = os.path.splitext(output_path)[1].lower()

    if ext in ['.yaml', '.yml']:
        with open(output_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    elif ext == '.json':
        with open(output_path, 'w') as f:
            json.dump(config, f, indent=2)
    else:
        raise ValueError(f"Unsupported config format: {ext}. Use .yaml, .yml, or .json")


def create_default_config() -> Dict[str, Any]:
    """
    Create default pipeline configuration.

    Returns:
        Default configuration dictionary
    """
    return {
        'antigen_pdb': 'path/to/antigen.pdb',
        'framework_pdb': 'path/to/framework.pdb',  # Optional
        'output_dir': './results',

        'design': {
            'mode': 'cdr_h3',
            'num_designs': 10,
            'num_steps': 50,
        },

        'epitope': {
            'residues': ['A10', 'A15', 'A20'],  # Optional
        },

        'docking': {
            'enabled': True,
            'hdock_bin': './bin/hdock',
            'createpl_bin': './bin/createpl',
        },

        'rfdiffusion': {
            'path': './RFdiffusion',
            'weights': None,  # Optional
        },

        'iglm': {
            'enabled': True,
            'species': '[HUMAN]',
            'temperature': 1.0,
        },

        'filtering': {
            'max_clashes': 100,
            'min_plddt': 70.0,
            'min_pdockq': 0.23,
            'max_sap_score': 50.0,
            'min_cdr_interface_pct': 60.0,
        },

        'structure_prediction': {
            'predictor': 'chai',
            'chai_path': None,
            'af3_path': None,
        },

        'device': 'cuda:0',
    }