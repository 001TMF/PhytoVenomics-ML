#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ConfigManager Class

This class handles loading, parsing, and managing configuration settings
for the Phytovenomics setup process.
"""

import os
import yaml
from pathlib import Path


class ConfigManager:
    """
    Handles configuration settings for the Phytovenomics setup process.
    
    Attributes:
        config_path (str): Path to the configuration file
        environment (str): Target environment (development, testing, production)
        config (dict): Loaded configuration data
    """
    
    def __init__(self, config_path='config/setup_config.yaml', environment='development'):
        """
        Initialize the ConfigManager.
        
        Args:
            config_path (str): Path to configuration file
            environment (str): Target environment (development, testing, production)
        """
        self.config_path = config_path
        self.environment = environment
        self.config = {}
        
        # Create config directory if it doesn't exist
        Path(os.path.dirname(config_path)).mkdir(exist_ok=True, parents=True)
    
    def load_config(self):
        """
        Load configuration from file.
        
        If the configuration file doesn't exist, generate a default configuration.
        
        Returns:
            dict: Configuration data
        """
        try:
            if not os.path.exists(self.config_path):
                self.generate_default_config()
            
            with open(self.config_path, 'r') as f:
                self.config = yaml.safe_load(f) or {}
            
            return self.config
        except Exception as e:
            print(f"Error loading configuration: {str(e)}")
            # Fall back to default configuration
            self.config = self._get_default_config()
            return self.config
    
    def get_environment_config(self):
        """
        Get configuration for the current environment.
        
        Returns:
            dict: Environment-specific configuration
        """
        if not self.config:
            self.load_config()
        
        # Get environment-specific configuration
        env_config = self.config.get('environments', {}).get(self.environment, {})
        
        # Merge with common configuration
        common_config = self.config.get('common', {})
        merged_config = {**common_config, **env_config}
        
        return merged_config
    
    def get_model_config(self):
        """
        Get model-related configuration.
        
        Returns:
            dict: Model configuration
        """
        if not self.config:
            self.load_config()
        
        return self.config.get('models', {})
    
    def update_config(self, key, value):
        """
        Update a configuration value.
        
        Args:
            key (str): Configuration key
            value: New value
            
        Returns:
            bool: Success status
        """
        if not self.config:
            self.load_config()
        
        # Support nested keys with dot notation
        if '.' in key:
            parts = key.split('.')
            curr = self.config
            for part in parts[:-1]:
                if part not in curr:
                    curr[part] = {}
                curr = curr[part]
            curr[parts[-1]] = value
        else:
            self.config[key] = value
        
        return True
    
    def save_config(self):
        """
        Save configuration to file.
        
        Returns:
            bool: Success status
        """
        try:
            with open(self.config_path, 'w') as f:
                yaml.dump(self.config, f, default_flow_style=False)
            return True
        except Exception as e:
            print(f"Error saving configuration: {str(e)}")
            return False
    
    def generate_default_config(self):
        """
        Generate default configuration file.
        
        Returns:
            bool: Success status
        """
        try:
            # Get default configuration
            default_config = self._get_default_config()
            
            # Save to file
            with open(self.config_path, 'w') as f:
                yaml.dump(default_config, f, default_flow_style=False)
            
            self.config = default_config
            return True
        except Exception as e:
            print(f"Error generating default configuration: {str(e)}")
            return False
    
    def _get_default_config(self):
        """
        Get default configuration.
        
        Returns:
            dict: Default configuration
        """
        return {
            'version': '1.0.0',
            'common': {
                'python_version': '>=3.8.0',
                'setup_date': None,
                'setup_completed': False
            },
            'environments': {
                'development': {
                    'dependencies': {
                        'core': 'requirements.txt',
                        'dev': 'requirements-dev.txt'
                    },
                    'required_models': [
                        'esm2_t33_650M_UR50D',
                        'igfold_base'
                    ],
                    'gpu_required': False,
                    'min_memory': '8GB',
                    'min_disk_space': '10GB'
                },
                'testing': {
                    'dependencies': {
                        'core': 'requirements.txt',
                        'test': 'requirements-test.txt'
                    },
                    'required_models': [
                        'esm2_t33_650M_UR50D',
                        'igfold_base'
                    ],
                    'gpu_required': False,
                    'min_memory': '8GB',
                    'min_disk_space': '10GB'
                },
                'production': {
                    'dependencies': {
                        'core': 'requirements.txt'
                    },
                    'required_models': [
                        'esm2_t36_3B_UR50D',
                        'igfold_base',
                        'igfold_improved'
                    ],
                    'gpu_required': True,
                    'min_memory': '16GB',
                    'min_disk_space': '50GB'
                }
            },
            'models': {
                'download_base_url': 'https://phytovenomics-models.s3.amazonaws.com',
                'default_download_dir': 'models',
                'available_models': {
                    'esm2_t33_650M_UR50D': {
                        'size': '1.3GB',
                        'version': '1.0',
                        'description': 'ESM-2 language model (650M parameters)',
                        'url': 'esm2_t33_650M_UR50D.pt',
                        'md5': 'e3846f627c3bc735ef8f6b2b89357240',
                        'required_disk_space': '2GB'
                    },
                    'esm2_t36_3B_UR50D': {
                        'size': '5.8GB',
                        'version': '1.0',
                        'description': 'ESM-2 language model (3B parameters)',
                        'url': 'esm2_t36_3B_UR50D.pt',
                        'md5': '9b2a4a6f0c7f45b4c51b1bf1d3cc66f7',
                        'required_disk_space': '12GB'
                    },
                    'igfold_base': {
                        'size': '450MB',
                        'version': '1.0',
                        'description': 'IgFold base model for antibody structure prediction',
                        'url': 'igfold_base.pt',
                        'md5': '134a32fb0e5d08a43ce89e09e54c680b',
                        'required_disk_space': '1GB'
                    },
                    'igfold_improved': {
                        'size': '850MB',
                        'version': '1.2',
                        'description': 'IgFold improved model with higher accuracy',
                        'url': 'igfold_improved.pt',
                        'md5': '9d15e03b0b78de8f880b24d60e8863e2',
                        'required_disk_space': '2GB'
                    }
                }
            },
            'model_download_dir': 'models',
            'log_dir': 'logs'
        }
