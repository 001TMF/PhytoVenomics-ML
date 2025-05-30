#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ModelManager Class

This class handles downloading, verification, and management of machine learning models
used by the Phytovenomics platform.
"""

import os
import sys
import hashlib
import requests
import shutil
import tempfile
from pathlib import Path
from tqdm import tqdm


class ModelManager:
    """
    Handles downloading and management of machine learning models.
    
    Attributes:
        model_config (dict): Configuration for models
        download_dir (str): Directory to download models to
        downloaded_models (list): List of successfully downloaded models
    """
    
    def __init__(self, model_config=None, download_dir='models'):
        """
        Initialize the ModelManager.
        
        Args:
            model_config (dict): Configuration for models
            download_dir (str): Directory to download models to
        """
        self.model_config = model_config or {}
        self.download_dir = download_dir
        self.downloaded_models = []
        
        # Create download directory if it doesn't exist
        Path(download_dir).mkdir(exist_ok=True, parents=True)
    
    def download_models(self, models_list=None, force=False):
        """
        Download specified models or all available models.
        
        Args:
            models_list (list): List of models to download, or None for all models
            force (bool): Whether to force redownload of existing models
            
        Returns:
            bool: Success status
        """
        # Get available models
        available_models = self.get_available_models()
        
        # If no models specified, download all available models
        if not models_list:
            models_list = available_models
        
        # Reset downloaded models list
        self.downloaded_models = []
        
        # Track success status
        success = True
        
        # Download each model
        for model_name in models_list:
            if model_name not in available_models:
                print(f"Warning: Model '{model_name}' not found in available models. Skipping.")
                continue
            
            # Download model
            result = self.download_model(model_name, force=force)
            if result:
                self.downloaded_models.append(model_name)
            else:
                success = False
        
        return success
    
    def download_model(self, model_name, force=False):
        """
        Download a specific model.
        
        Args:
            model_name (str): Name of the model to download
            force (bool): Whether to force redownload if the model already exists
            
        Returns:
            bool: Success status
        """
        # Check if model exists in configuration
        if model_name not in self.model_config.get('available_models', {}):
            print(f"Error: Model '{model_name}' not found in available models.")
            return False
        
        # Get model information
        model_info = self.model_config['available_models'][model_name]
        model_url = f"{self.model_config.get('download_base_url', '')}/{model_info['url']}"
        model_path = os.path.join(self.download_dir, model_info['url'])
        
        # Check if model already exists and is valid
        if os.path.exists(model_path) and not force:
            if self.verify_model_checksum(model_name):
                print(f"Model '{model_name}' already exists and is valid. Skipping.")
                return True
            else:
                print(f"Model '{model_name}' exists but is invalid. Redownloading.")
        
        # Create temporary directory for download
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = os.path.join(temp_dir, model_info['url'])
            
            # Download the model
            print(f"Downloading model '{model_name}' ({model_info['size']})...")
            try:
                # Send a HEAD request to get the file size
                response = requests.head(model_url)
                total_size = int(response.headers.get('content-length', 0))
                
                # Download with progress bar
                with requests.get(model_url, stream=True) as r:
                    r.raise_for_status()
                    with open(temp_path, 'wb') as f, tqdm(
                        desc=model_name,
                        total=total_size,
                        unit='B',
                        unit_scale=True,
                        unit_divisor=1024,
                    ) as pbar:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                pbar.update(len(chunk))
                
                # Create model directory if needed
                os.makedirs(os.path.dirname(model_path), exist_ok=True)
                
                # Move file to final location
                shutil.move(temp_path, model_path)
                
                # Verify download
                if self.verify_model_checksum(model_name):
                    print(f"Successfully downloaded and verified model '{model_name}'.")
                    
                    # Extract model if needed
                    if model_path.endswith('.tar.gz') or model_path.endswith('.zip'):
                        self.extract_model(model_path)
                    
                    return True
                else:
                    print(f"Error: Failed to verify model '{model_name}'. Checksum mismatch.")
                    if os.path.exists(model_path):
                        os.remove(model_path)
                    return False
                
            except Exception as e:
                print(f"Error downloading model '{model_name}': {str(e)}")
                if os.path.exists(model_path):
                    os.remove(model_path)
                return False
    
    def verify_model_checksum(self, model_name):
        """
        Verify the checksum of a downloaded model.
        
        Args:
            model_name (str): Name of the model to verify
            
        Returns:
            bool: Whether the checksum matches
        """
        # Check if model exists in configuration
        if model_name not in self.model_config.get('available_models', {}):
            print(f"Error: Model '{model_name}' not found in available models.")
            return False
        
        # Get model information
        model_info = self.model_config['available_models'][model_name]
        expected_md5 = model_info.get('md5', '')
        if not expected_md5:
            print(f"Warning: No MD5 checksum provided for model '{model_name}'. Skipping verification.")
            return True
        
        # Get model path
        model_path = os.path.join(self.download_dir, model_info['url'])
        if not os.path.exists(model_path):
            print(f"Error: Model file not found at {model_path}")
            return False
        
        # Calculate MD5 checksum
        md5_hash = hashlib.md5()
        with open(model_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b""):
                md5_hash.update(chunk)
        calculated_md5 = md5_hash.hexdigest()
        
        # Compare checksums
        if calculated_md5 == expected_md5:
            return True
        else:
            print(f"Checksum mismatch for model '{model_name}'")
            print(f"Expected: {expected_md5}")
            print(f"Actual: {calculated_md5}")
            return False
    
    def extract_model(self, model_path):
        """
        Extract a compressed model file.
        
        Args:
            model_path (str): Path to the model file
            
        Returns:
            bool: Success status
        """
        try:
            extract_dir = os.path.dirname(model_path)
            print(f"Extracting {model_path} to {extract_dir}...")
            
            if model_path.endswith('.tar.gz'):
                import tarfile
                with tarfile.open(model_path, 'r:gz') as tar:
                    tar.extractall(path=extract_dir)
                return True
            elif model_path.endswith('.zip'):
                import zipfile
                with zipfile.ZipFile(model_path, 'r') as zip_ref:
                    zip_ref.extractall(extract_dir)
                return True
            else:
                print(f"Unsupported archive format for {model_path}")
                return False
        except Exception as e:
            print(f"Error extracting model: {str(e)}")
            return False
    
    def get_model_size(self, model_name):
        """
        Get the size of a model.
        
        Args:
            model_name (str): Name of the model
            
        Returns:
            int: Size in bytes, or 0 if unknown
        """
        # Check if model exists in configuration
        if model_name not in self.model_config.get('available_models', {}):
            print(f"Error: Model '{model_name}' not found in available models.")
            return 0
        
        # Get model information
        model_info = self.model_config['available_models'][model_name]
        size_str = model_info.get('size', '0')
        
        # Parse size string (e.g., "1.3GB")
        try:
            size_num = float(size_str[:-2])
            size_unit = size_str[-2:]
            
            if size_unit == 'KB':
                return int(size_num * 1024)
            elif size_unit == 'MB':
                return int(size_num * 1024 * 1024)
            elif size_unit == 'GB':
                return int(size_num * 1024 * 1024 * 1024)
            else:
                return 0
        except:
            return 0
    
    def get_available_models(self):
        """
        Get a list of available models.
        
        Returns:
            list: List of available model names
        """
        return list(self.model_config.get('available_models', {}).keys())
    
    def get_downloaded_models(self):
        """
        Get a list of successfully downloaded models.
        
        Returns:
            list: List of downloaded model names
        """
        return self.downloaded_models
