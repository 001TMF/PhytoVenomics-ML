#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
SetupManager Class

This class orchestrates the setup process for the Phytovenomics project,
managing environment setup, dependency installation, model downloading,
and project configuration.
"""

import os
import sys
import time
import logging
from pathlib import Path

from .config_manager import ConfigManager
from .environment_manager import EnvironmentManager
from .dependency_manager import DependencyManager
from .model_manager import ModelManager
from .system_validator import SystemValidator
from .logger import Logger


class SetupManager:
    """
    Orchestrates the setup process for the Phytovenomics project.
    
    Attributes:
        config_path (str): Path to the configuration file
        environment (str): Target environment (development, testing, production)
        force (bool): Whether to force reinstallation/redownload
        logger (Logger): Logger instance
    """
    
    def __init__(self, config_path='config/setup_config.yaml', 
                 environment='development', force=False):
        """
        Initialize the SetupManager.
        
        Args:
            config_path (str): Path to configuration file
            environment (str): Target environment (development, testing, production)
            force (bool): Whether to force reinstallation/redownload
        """
        self.config_path = config_path
        self.environment = environment
        self.force = force
        
        # Initialize logger
        log_dir = Path('logs')
        log_dir.mkdir(exist_ok=True)
        log_path = log_dir / f"setup_{time.strftime('%Y%m%d_%H%M%S')}.log"
        self.logger = Logger(log_path=str(log_path), log_level="INFO")
        
        # Load configuration
        self.config_manager = ConfigManager(config_path=config_path, 
                                           environment=environment)
        self.config = self.config_manager.load_config()
        
        # Initialize component managers
        self.system_validator = SystemValidator()
        self.env_manager = EnvironmentManager(environment=environment, 
                                             config=self.config)
        self.dependency_manager = DependencyManager(environment=environment, 
                                                  config=self.config)
        self.model_manager = ModelManager(
            model_config=self.config_manager.get_model_config(),
            download_dir=self.config.get('model_download_dir', 'models')
        )
        
        # Setup report tracking
        self.setup_report = {
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'environment': environment,
            'system_validation': {},
            'environment_setup': False,
            'dependencies_installed': [],
            'models_downloaded': [],
            'tests_passed': False
        }
    
    def run_setup(self, skip_tests=False):
        """
        Run the complete setup process.
        
        Args:
            skip_tests (bool): Whether to skip running tests after setup
            
        Returns:
            bool: Success status
        """
        self.logger.log_info(f"Starting Phytovenomics setup for {self.environment} environment")
        
        # Step 1: Validate system requirements
        if not self.validate_system_requirements():
            self.logger.log_error("System requirements not met. Setup failed.")
            return False
        
        # Step 2: Setup environment
        if not self.setup_environment():
            self.logger.log_error("Environment setup failed.")
            return False
        
        # Step 3: Install dependencies
        if not self.install_dependencies():
            self.logger.log_error("Dependency installation failed.")
            return False
        
        # Step 4: Download models
        if not self.download_models():
            self.logger.log_error("Model download failed.")
            return False
        
        # Step 5: Configure project
        if not self.configure_project():
            self.logger.log_error("Project configuration failed.")
            return False
        
        # Step 6: Run tests (optional)
        if not skip_tests:
            if not self.run_tests():
                self.logger.log_warning("Some tests failed.")
        else:
            self.logger.log_info("Skipping tests as requested.")
        
        # Update config to reflect successful setup
        self.config_manager.update_config("setup_completed", True)
        self.config_manager.update_config("last_setup_time", time.strftime('%Y-%m-%d %H:%M:%S'))
        self.config_manager.save_config()
        
        self.logger.log_success("Phytovenomics setup completed successfully!")
        return True
    
    def validate_system_requirements(self):
        """
        Validate system requirements for the project.
        
        Returns:
            bool: Whether the system meets all requirements
        """
        self.logger.log_info("Validating system requirements...")
        
        # Run system validation
        validation_results = self.system_validator.validate_system()
        self.setup_report['system_validation'] = validation_results
        
        if all(validation_results.values()):
            self.logger.log_success("System meets all requirements.")
            return True
        else:
            self.logger.log_error("System does not meet requirements:")
            for check, result in validation_results.items():
                if not result:
                    self.logger.log_error(f"- {check}: Failed")
            return False
    
    def setup_environment(self):
        """
        Set up the project environment.
        
        Returns:
            bool: Success status
        """
        self.logger.log_info("Setting up environment...")
        
        # Create and activate environment
        result = self.env_manager.setup_environment()
        self.setup_report['environment_setup'] = result
        
        if result:
            self.logger.log_success("Environment setup completed.")
        else:
            self.logger.log_error("Environment setup failed.")
        
        return result
    
    def install_dependencies(self):
        """
        Install project dependencies.
        
        Returns:
            bool: Success status
        """
        self.logger.log_info("Installing dependencies...")
        
        # Install core dependencies
        core_result = self.dependency_manager.install_core_dependencies()
        if not core_result:
            self.logger.log_error("Failed to install core dependencies.")
            return False
        
        # Install environment-specific dependencies
        env_result = self.dependency_manager.install_environment_specific_dependencies()
        if not env_result:
            self.logger.log_error("Failed to install environment-specific dependencies.")
            return False
        
        # Generate dependency report
        dep_report = self.dependency_manager.generate_dependency_report()
        self.setup_report['dependencies_installed'] = dep_report['installed_packages']
        
        self.logger.log_success("Dependency installation completed.")
        return True
    
    def download_models(self):
        """
        Download required models.
        
        Returns:
            bool: Success status
        """
        self.logger.log_info("Downloading required models...")
        
        # Get required models for this environment
        env_config = self.config_manager.get_environment_config()
        required_models = env_config.get('required_models', [])
        
        # Download models
        result = self.model_manager.download_models(models_list=required_models, 
                                                 force=self.force)
        
        if result:
            self.setup_report['models_downloaded'] = self.model_manager.get_downloaded_models()
            self.logger.log_success("Model download completed.")
        else:
            self.logger.log_error("Model download failed.")
        
        return result
    
    def configure_project(self):
        """
        Configure the project after installation.
        
        Returns:
            bool: Success status
        """
        self.logger.log_info("Configuring project...")
        
        try:
            # Create necessary directories
            for directory in ['results', 'logs', 'data/temp']:
                Path(directory).mkdir(exist_ok=True, parents=True)
            
            # Configure environment-specific settings
            env_config = self.config_manager.get_environment_config()
            
            # Create environment-specific configuration file
            env_config_path = f"config/env_{self.environment}.yaml"
            if not os.path.exists(env_config_path) or self.force:
                with open(env_config_path, 'w') as f:
                    import yaml
                    yaml.dump(env_config, f, default_flow_style=False)
            
            self.logger.log_success("Project configuration completed.")
            return True
        except Exception as e:
            self.logger.log_error(f"Project configuration failed: {str(e)}")
            return False
    
    def run_tests(self):
        """
        Run tests to verify the setup.
        
        Returns:
            bool: Whether all tests passed
        """
        self.logger.log_info("Running tests...")
        
        try:
            # Import test module conditionally to avoid circular imports
            from phytovenomics.setup.test_runner import TestRunner
            
            # Run tests
            test_runner = TestRunner()
            test_results = test_runner.run_setup_tests(environment=self.environment)
            
            # Update report
            self.setup_report['tests_passed'] = test_results['all_passed']
            
            if test_results['all_passed']:
                self.logger.log_success("All tests passed.")
                return True
            else:
                self.logger.log_warning(f"Some tests failed: {test_results['failed_tests']}")
                return False
        except Exception as e:
            self.logger.log_error(f"Error running tests: {str(e)}")
            return False
    
    def generate_setup_report(self):
        """
        Generate a report of the setup process.
        
        Returns:
            dict: Setup report
        """
        return self.setup_report
