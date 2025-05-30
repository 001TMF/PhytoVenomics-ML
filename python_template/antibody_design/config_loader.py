#!/usr/bin/env python
# antibody_design/config_loader.py

import os
import yaml
import logging
from typing import Dict, Any, Union
from pathlib import Path

logger = logging.getLogger(__name__)

def load_config(config_path_or_dict: Union[str, Dict, Path]) -> Dict[str, Any]:
    """
    Load configuration from a YAML file or dictionary.
    
    Args:
        config_path_or_dict: Path to YAML configuration file or dictionary with configuration
        
    Returns:
        Dictionary with configuration parameters
    
    Raises:
        FileNotFoundError: If the configuration file doesn't exist
        ValueError: If the configuration file is malformed
    """
    # If input is already a dictionary, return it
    if isinstance(config_path_or_dict, dict):
        logger.debug("Configuration provided as dictionary")
        return config_path_or_dict
    
    # Convert to Path object if it's a string
    if isinstance(config_path_or_dict, str):
        config_path = Path(config_path_or_dict)
    else:
        config_path = config_path_or_dict
    
    # Check if file exists
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found at {config_path}")
    
    # Load YAML file
    try:
        logger.info(f"Loading configuration from {config_path}")
        with open(config_path, 'r') as config_file:
            config = yaml.safe_load(config_file)
        
        if not config:
            raise ValueError(f"Empty or invalid configuration file: {config_path}")
            
        # Extract nested sections to top level if needed
        if "general" in config:
            for key, value in config["general"].items():
                # Only set if not already in the top level
                if key not in config:
                    config[key] = value
                    
        # Extract evolution parameters
        if "evolution" in config:
            for key, value in config["evolution"].items():
                if key not in config:
                    config[key] = value
        
        # Extract objectives and weights
        if "objectives" in config:
            config["objective_weights"] = config["objectives"]
            
            # Also ensure there's an objectives list for algorithms that expect it
            if "objective_weights" in config and isinstance(config["objective_weights"], dict):
                config["objectives"] = list(config["objective_weights"].keys())
            
        return config
        
    except yaml.YAMLError as e:
        logger.error(f"Error parsing YAML configuration: {e}")
        raise ValueError(f"Invalid YAML in configuration file: {e}")
    except Exception as e:
        logger.error(f"Error loading configuration: {e}")
        raise