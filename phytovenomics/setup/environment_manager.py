#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
EnvironmentManager Class

This class handles creating and configuring the Python environment
for the Phytovenomics platform, including virtual environment setup.
"""

import os
import sys
import subprocess
import platform
from pathlib import Path


class EnvironmentManager:
    """
    Handles setting up and configuring the Python environment.
    
    Attributes:
        environment (str): Target environment (development, testing, production)
        config (dict): Configuration data
    """
    
    def __init__(self, environment='development', config=None):
        """
        Initialize the EnvironmentManager.
        
        Args:
            environment (str): Target environment (development, testing, production)
            config (dict): Configuration data
        """
        self.environment = environment
        self.config = config or {}
    
    def setup_environment(self):
        """
        Set up the Python environment.
        
        Returns:
            bool: Success status
        """
        print(f"Setting up {self.environment} environment...")
        
        # Get environment-specific configuration
        env_config = self.config.get('environments', {}).get(self.environment, {})
        
        # Create virtual environment if not already in one
        if not self._is_in_virtualenv():
            venv_name = f"phytovenomics_{self.environment}"
            success = self.create_virtual_env(venv_name)
            if not success:
                print("Failed to create virtual environment.")
                return False
            
            # Activate environment
            success = self.activate_environment(venv_name)
            if not success:
                print("Failed to activate virtual environment.")
                return False
        
        # Verify environment
        return self.verify_environment()
    
    def create_virtual_env(self, name):
        """
        Create a virtual environment.
        
        Args:
            name (str): Name of the virtual environment
            
        Returns:
            bool: Success status
        """
        try:
            print(f"Creating virtual environment '{name}'...")
            
            # Create virtual environment directory if it doesn't exist
            venv_dir = Path('.venv') / name
            if not venv_dir.exists():
                venv_dir.parent.mkdir(exist_ok=True)
                
                # Use the appropriate command for creating a virtual environment
                import venv
                venv.create(venv_dir, with_pip=True)
                
                print(f"Virtual environment created at {venv_dir}")
            else:
                print(f"Virtual environment already exists at {venv_dir}")
            
            return True
        except Exception as e:
            print(f"Error creating virtual environment: {str(e)}")
            return False
    
    def activate_environment(self, venv_name=None):
        """
        Activate the virtual environment.
        
        Note: In Python scripts, we don't actually "activate" the environment in the
        traditional sense, but we print instructions for the user to do so.
        
        Args:
            venv_name (str): Name of the virtual environment
            
        Returns:
            bool: Success status (always True, since we're just printing instructions)
        """
        venv_dir = Path('.venv') / (venv_name or f"phytovenomics_{self.environment}")
        
        # Print activation instructions based on the platform
        if venv_dir.exists():
            print("\nTo activate the virtual environment, run:")
            if platform.system() == "Windows":
                print(f"> {venv_dir}\\Scripts\\activate")
            else:
                print(f"$ source {venv_dir}/bin/activate")
        
        return True
    
    def verify_environment(self):
        """
        Verify that the environment is properly set up.
        
        Returns:
            bool: Whether the environment is properly set up
        """
        try:
            # Check Python version
            import platform
            python_version = platform.python_version()
            print(f"Python version: {python_version}")
            
            # Check if pip is available
            subprocess.check_call([sys.executable, '-m', 'pip', '--version'],
                                stdout=subprocess.DEVNULL, 
                                stderr=subprocess.DEVNULL)
            
            # Check if we're in a virtual environment
            in_venv = self._is_in_virtualenv()
            if in_venv:
                print(f"Running in virtual environment: {sys.prefix}")
            else:
                print("Warning: Not running in a virtual environment.")
            
            return True
        except Exception as e:
            print(f"Error verifying environment: {str(e)}")
            return False
    
    def _is_in_virtualenv(self):
        """
        Check if we're currently running in a virtual environment.
        
        Returns:
            bool: Whether we're in a virtual environment
        """
        # Check if we're in a virtual environment
        return (hasattr(sys, 'real_prefix') or 
                (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix))
