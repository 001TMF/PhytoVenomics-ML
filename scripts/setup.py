#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Phytovenomics Setup Script

This script initializes the Phytovenomics environment, installs dependencies,
downloads required models, and configures the project for use.

Usage:
    python setup.py [OPTIONS]

Options:
    --environment      The environment to set up (development, testing, production)
    --config           Path to configuration file
    --skip-tests       Skip running tests after setup
    --force            Force reinstallation of dependencies and models
    --models           Specific models to download (comma-separated)
    --help             Show this help message and exit
"""

import os
import sys
import click
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from phytovenomics.setup.setup_manager import SetupManager
from phytovenomics.setup.model_manager import ModelManager
from phytovenomics.setup.config_manager import ConfigManager
from phytovenomics.setup.system_validator import SystemValidator


@click.group()
def cli():
    """Phytovenomics setup commands."""
    pass


@cli.command()
@click.option('--environment', '-e', default='development',
              type=click.Choice(['development', 'testing', 'production']),
              help='Environment to setup (development, testing, production)')
@click.option('--config', '-c', default='config/setup_config.yaml',
              help='Path to configuration file')
@click.option('--skip-tests', is_flag=True, help='Skip running tests after setup')
@click.option('--force', '-f', is_flag=True, help='Force reinstallation of dependencies and models')
def setup(environment, config, skip_tests, force):
    """Set up the Phytovenomics environment."""
    click.echo(f"Setting up Phytovenomics for {environment} environment")
    
    # Initialize setup manager
    setup_manager = SetupManager(
        config_path=config,
        environment=environment,
        force=force
    )
    
    # Run setup
    result = setup_manager.run_setup(skip_tests=skip_tests)
    
    if result:
        click.echo(click.style("✅ Phytovenomics setup completed successfully!", fg='green'))
        
        # Show setup report
        report = setup_manager.generate_setup_report()
        click.echo("\nSetup Summary:")
        click.echo(f"Environment: {environment}")
        click.echo(f"Dependencies installed: {report['dependencies_installed']}")
        click.echo(f"Models downloaded: {', '.join(report['models_downloaded'])}")
        
        click.echo("\nNext steps:")
        click.echo("  1. Run 'python demo_test.py' to verify the installation")
        click.echo("  2. Review the documentation at docs/")
    else:
        click.echo(click.style("❌ Phytovenomics setup failed.", fg='red'))
        click.echo("Check the logs for more information.")
        sys.exit(1)


@cli.command()
@click.option('--models', '-m', help='Comma-separated list of models to download')
@click.option('--force', '-f', is_flag=True, help='Force redownload of models')
@click.option('--config', '-c', default='config/setup_config.yaml',
              help='Path to configuration file')
def download_models(models, force, config):
    """Download specific ML models."""
    # Load config
    config_manager = ConfigManager(config_path=config)
    model_config = config_manager.get_model_config()
    
    # Initialize model manager
    model_manager = ModelManager(model_config=model_config)
    
    # Parse models list
    models_list = None
    if models:
        models_list = [m.strip() for m in models.split(',')]
        click.echo(f"Downloading specified models: {', '.join(models_list)}")
    else:
        click.echo("Downloading all required models")
    
    # Download models
    result = model_manager.download_models(models_list=models_list, force=force)
    
    if result:
        click.echo(click.style("✅ Model download completed successfully!", fg='green'))
    else:
        click.echo(click.style("❌ Model download failed.", fg='red'))
        click.echo("Check the logs for more information.")
        sys.exit(1)


@cli.command()
def check_system():
    """Check if the system meets the requirements."""
    click.echo("Checking system requirements...")
    
    # Initialize system validator
    validator = SystemValidator()
    
    # Run validation
    results = validator.validate_system(verbose=True)
    
    if all(results.values()):
        click.echo(click.style("✅ System meets all requirements!", fg='green'))
    else:
        click.echo(click.style("❌ System does not meet all requirements:", fg='yellow'))
        for check, passed in results.items():
            status = "✓" if passed else "✗"
            color = "green" if passed else "red"
            click.echo(click.style(f"{status} {check}", fg=color))


@cli.command()
@click.option('--environment', '-e', default='development',
              type=click.Choice(['development', 'testing', 'production']),
              help='Environment to update dependencies for')
@click.option('--config', '-c', default='config/setup_config.yaml',
              help='Path to configuration file')
def update_dependencies(environment, config):
    """Update project dependencies."""
    click.echo(f"Updating dependencies for {environment} environment")
    
    # Initialize setup manager
    setup_manager = SetupManager(
        config_path=config,
        environment=environment,
        force=True
    )
    
    # Update dependencies
    result = setup_manager.install_dependencies()
    
    if result:
        click.echo(click.style("✅ Dependencies updated successfully!", fg='green'))
    else:
        click.echo(click.style("❌ Dependency update failed.", fg='red'))
        click.echo("Check the logs for more information.")
        sys.exit(1)


if __name__ == "__main__":
    cli()
