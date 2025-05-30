#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DependencyManager Class

This class handles installing and managing dependencies for the Phytovenomics platform,
including core, environment-specific, and optional dependencies.
"""

import os
import sys
import subprocess
import pkg_resources
from pathlib import Path


class DependencyManager:
    """
    Handles installing and managing project dependencies.
    
    Attributes:
        environment (str): Target environment (development, testing, production)
        config (dict): Configuration data
        requirements_files (dict): Dictionary of requirements file paths
    """
    
    def __init__(self, environment='development', config=None):
        """
        Initialize the DependencyManager.
        
        Args:
            environment (str): Target environment (development, testing, production)
            config (dict): Configuration data
        """
        self.environment = environment
        self.config = config or {}
        
        # Get environment-specific configuration
        env_config = self.config.get('environments', {}).get(self.environment, {})
        self.requirements_files = env_config.get('dependencies', {})
    
    def install_dependencies(self):
        """
        Install all required dependencies.
        
        Returns:
            bool: Success status
        """
        # Install core dependencies
        if not self.install_core_dependencies():
            return False
        
        # Install environment-specific dependencies
        if not self.install_environment_specific_dependencies():
            return False
        
        # Verify installations
        if not self.verify_installations():
            return False
        
        return True
    
    def install_core_dependencies(self):
        """
        Install core dependencies from requirements file.
        
        Returns:
            bool: Success status
        """
        # Get core requirements file path
        req_file = self.requirements_files.get('core')
        if not req_file:
            print("Warning: No core requirements file specified.")
            return True
        
        # Install requirements
        return self._install_from_requirements(req_file)
    
    def install_environment_specific_dependencies(self):
        """
        Install environment-specific dependencies.
        
        Returns:
            bool: Success status
        """
        # Get environment-specific requirements
        env_specific_req = None
        if self.environment == 'development':
            env_specific_req = self.requirements_files.get('dev')
        elif self.environment == 'testing':
            env_specific_req = self.requirements_files.get('test')
        elif self.environment == 'production':
            env_specific_req = self.requirements_files.get('prod')
        
        if not env_specific_req:
            print(f"No specific requirements file for {self.environment} environment.")
            return True
        
        # Install requirements
        return self._install_from_requirements(env_specific_req)
    
    def install_optional_dependencies(self, category=None):
        """
        Install optional dependencies.
        
        Args:
            category (str): Category of optional dependencies to install
            
        Returns:
            bool: Success status
        """
        # Get optional requirements
        optional_req = self.requirements_files.get('optional', {})
        if not optional_req:
            print("No optional requirements specified.")
            return True
        
        if category:
            # Install specific category
            req_file = optional_req.get(category)
            if not req_file:
                print(f"No requirements file for optional category '{category}'.")
                return False
            
            return self._install_from_requirements(req_file)
        else:
            # Install all optional categories
            all_success = True
            for category, req_file in optional_req.items():
                success = self._install_from_requirements(req_file)
                all_success = all_success and success
            
            return all_success
    
    def verify_installations(self):
        """
        Verify that all required packages are installed.
        
        Returns:
            bool: Whether all required packages are installed
        """
        try:
            # Get all requirements files
            req_files = []
            for key, value in self.requirements_files.items():
                if isinstance(value, str):
                    req_files.append(value)
                elif isinstance(value, dict):
                    req_files.extend(value.values())
            
            # Check each requirements file
            all_installed = True
            for req_file in req_files:
                if not self._verify_requirements_installed(req_file):
                    all_installed = False
            
            return all_installed
        except Exception as e:
            print(f"Error verifying installations: {str(e)}")
            return False
    
    def generate_dependency_report(self):
        """
        Generate a report of installed dependencies.
        
        Returns:
            dict: Report of installed packages
        """
        report = {
            'installed_packages': [],
            'missing_packages': []
        }
        
        # Get installed packages
        installed_packages = {pkg.key: pkg.version for pkg in pkg_resources.working_set}
        report['installed_packages'] = [f"{key}=={version}" for key, version in installed_packages.items()]
        
        # Get all requirements files
        req_files = []
        for key, value in self.requirements_files.items():
            if isinstance(value, str):
                req_files.append(value)
            elif isinstance(value, dict):
                req_files.extend(value.values())
        
        # Check for missing packages
        for req_file in req_files:
            try:
                with open(req_file, 'r') as f:
                    for line in f:
                        # Parse requirement line
                        line = line.strip()
                        if not line or line.startswith('#'):
                            continue
                        
                        # Split package name and version
                        if '==' in line:
                            pkg_name = line.split('==')[0].lower()
                        elif '>=' in line:
                            pkg_name = line.split('>=')[0].lower()
                        elif '<=' in line:
                            pkg_name = line.split('<=')[0].lower()
                        else:
                            pkg_name = line.lower()
                        
                        # Check if package is installed
                        if pkg_name not in installed_packages:
                            report['missing_packages'].append(line)
            except Exception:
                pass
        
        return report
    
    def _install_from_requirements(self, req_file):
        """
        Install packages from a requirements file.
        
        Args:
            req_file (str): Path to requirements file
            
        Returns:
            bool: Success status
        """
        try:
            if not os.path.exists(req_file):
                print(f"Requirements file not found: {req_file}")
                return False
            
            print(f"Installing dependencies from {req_file}...")
            
            # Run pip install
            result = subprocess.run(
                [sys.executable, '-m', 'pip', 'install', '-r', req_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False
            )
            
            if result.returncode != 0:
                print(f"Error installing dependencies from {req_file}:")
                print(result.stderr)
                return False
            else:
                print(f"Successfully installed dependencies from {req_file}")
                return True
        except Exception as e:
            print(f"Error installing dependencies: {str(e)}")
            return False
    
    def _verify_requirements_installed(self, req_file):
        """
        Verify that all packages in a requirements file are installed.
        
        Args:
            req_file (str): Path to requirements file
            
        Returns:
            bool: Whether all packages are installed
        """
        try:
            if not os.path.exists(req_file):
                print(f"Requirements file not found: {req_file}")
                return False
            
            # Get installed packages
            installed_packages = {pkg.key: pkg.version for pkg in pkg_resources.working_set}
            
            # Check each requirement
            with open(req_file, 'r') as f:
                for line in f:
                    # Parse requirement line
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # Split package name and version
                    if '==' in line:
                        pkg_name = line.split('==')[0].lower()
                    elif '>=' in line:
                        pkg_name = line.split('>=')[0].lower()
                    elif '<=' in line:
                        pkg_name = line.split('<=')[0].lower()
                    else:
                        pkg_name = line.lower()
                    
                    # Check if package is installed
                    if pkg_name not in installed_packages:
                        print(f"Missing package: {pkg_name}")
                        return False
            
            return True
        except Exception as e:
            print(f"Error verifying requirements: {str(e)}")
            return False
