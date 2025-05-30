#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
SystemValidator Class

This class checks whether the system meets the requirements for the Phytovenomics platform,
including Python version, disk space, memory, and GPU availability.
"""

import os
import sys
import platform
import shutil
import psutil
from packaging import version


class SystemValidator:
    """
    Validates system requirements for the Phytovenomics platform.
    """
    
    def validate_system(self, verbose=False):
        """
        Validate all system requirements.
        
        Args:
            verbose (bool): Whether to print validation results
            
        Returns:
            dict: Validation results for each check
        """
        results = {
            'Python Version': self.check_python_version(),
            'Disk Space': self.check_disk_space(),
            'Memory': self.check_memory(),
            'GPU Available': self.check_gpu(),
            'CUDA Available': self.check_cuda(),
            'Internet Connection': self.check_internet_connection()
        }
        
        if verbose:
            print("System Validation Results:")
            for check, result in results.items():
                status = "✓" if result else "✗"
                print(f"{status} {check}")
        
        return results
    
    def check_python_version(self, min_version='3.8.0'):
        """
        Check if Python version meets requirements.
        
        Args:
            min_version (str): Minimum required Python version
            
        Returns:
            bool: Whether Python version meets requirements
        """
        current_version = platform.python_version()
        return version.parse(current_version) >= version.parse(min_version)
    
    def check_disk_space(self, min_gb=10):
        """
        Check if available disk space meets requirements.
        
        Args:
            min_gb (int): Minimum required disk space in GB
            
        Returns:
            bool: Whether disk space meets requirements
        """
        # Get disk usage of the current directory
        disk = shutil.disk_usage('.')
        free_gb = disk.free / (1024**3)  # Convert bytes to GB
        return free_gb >= min_gb
    
    def check_memory(self, min_gb=8):
        """
        Check if available memory meets requirements.
        
        Args:
            min_gb (int): Minimum required memory in GB
            
        Returns:
            bool: Whether memory meets requirements
        """
        mem_info = psutil.virtual_memory()
        total_gb = mem_info.total / (1024**3)  # Convert bytes to GB
        return total_gb >= min_gb
    
    def check_gpu(self):
        """
        Check if GPU is available.
        
        Returns:
            bool: Whether GPU is available
        """
        try:
            # Try to import torch and check if CUDA is available
            import torch
            return torch.cuda.is_available()
        except ImportError:
            # If torch is not installed, try to check for NVIDIA GPU using system commands
            try:
                if platform.system() == 'Windows':
                    # Use Windows Management Instrumentation (WMI)
                    import wmi
                    wmi_obj = wmi.WMI()
                    for gpu in wmi_obj.Win32_VideoController():
                        if 'nvidia' in gpu.Name.lower():
                            return True
                    return False
                elif platform.system() == 'Linux':
                    # Use nvidia-smi on Linux
                    return os.system('which nvidia-smi > /dev/null 2>&1') == 0
                else:
                    # macOS or other systems - less reliable check
                    return False
            except:
                return False
    
    def check_cuda(self, min_version='11.0'):
        """
        Check if CUDA is available and meets minimum version.
        
        Args:
            min_version (str): Minimum required CUDA version
            
        Returns:
            bool: Whether CUDA is available and meets minimum version
        """
        try:
            # Try to import torch and check CUDA version
            import torch
            if not torch.cuda.is_available():
                return False
            
            # Get CUDA version
            cuda_version = torch.version.cuda
            if cuda_version is None:
                return False
            
            return version.parse(cuda_version) >= version.parse(min_version)
        except ImportError:
            # If torch is not installed, try to check CUDA using system commands
            try:
                if platform.system() == 'Linux':
                    # Get CUDA version using nvcc
                    import subprocess
                    output = subprocess.check_output(['nvcc', '--version']).decode('utf-8')
                    for line in output.split('\n'):
                        if 'release' in line:
                            # Extract version number
                            import re
                            match = re.search(r'release (\d+\.\d+)', line)
                            if match:
                                cuda_version = match.group(1)
                                return version.parse(cuda_version) >= version.parse(min_version)
                return False
            except:
                return False
    
    def check_internet_connection(self):
        """
        Check if internet connection is available.
        
        Returns:
            bool: Whether internet connection is available
        """
        try:
            import urllib.request
            urllib.request.urlopen('http://google.com', timeout=1)
            return True
        except:
            return False
