# phytovenomics/cli/run_cli.py
#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Phytovenomics Run CLI

This module provides a command-line interface for running the Phytovenomics
analysis pipelines, models, and tools.
"""

import sys
import click
from pathlib import Path

from phytovenomics.setup.config_manager import ConfigManager


@click.group()
def cli():
    """Phytovenomics run commands."""
    pass


@cli.command()
@click.option('--config', '-c', default='config/setup_config.yaml',
              help='Path to configuration file')
def demo(config):
    """Run a demo to verify installation."""
    click.echo("Running Phytovenomics demonstration...")
    
    try:
        # Import the demo module
        from phytovenomics.demo import run_demo
        
        # Run the demo
        result = run_demo()
        
        if result:
            click.echo(click.style("✅ Demo completed successfully!", fg='green'))
            click.echo("Your Phytovenomics installation is working correctly.")
        else:
            click.echo(click.style("❌ Demo failed.", fg='red'))
            click.echo("Check the logs for more information.")
            sys.exit(1)
    except ImportError:
        click.echo(click.style("❌ Could not import demo module.", fg='red'))
        click.echo("Make sure Phytovenomics is properly installed.")
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"❌ Demo failed with error: {str(e)}", fg='red'))
        sys.exit(1)


@cli.command()
@click.argument('sequence_file', type=click.Path(exists=True))
@click.option('--model', '-m', default='default',
              help='Model to use for prediction')
@click.option('--output', '-o', default='results/predictions.json',
              help='Path to save prediction results')
@click.option('--config', '-c', default='config/setup_config.yaml',
              help='Path to configuration file')
def predict(sequence_file, model, output, config):
    """Run prediction on a sequence file."""
    click.echo(f"Running prediction on {sequence_file} using {model} model...")
    
    try:
        # Import the prediction module
        from phytovenomics.prediction import run_prediction
        
        # Ensure output directory exists
        output_dir = Path(output).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run prediction
        result = run_prediction(sequence_file, model=model, output_file=output)
        
        if result:
            click.echo(click.style(f"✅ Prediction completed successfully! Results saved to {output}", fg='green'))
        else:
            click.echo(click.style("❌ Prediction failed.", fg='red'))
            sys.exit(1)
    except ImportError:
        click.echo(click.style("❌ Could not import prediction module.", fg='red'))
        click.echo("Make sure Phytovenomics is properly installed.")
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"❌ Prediction failed with error: {str(e)}", fg='red'))
        sys.exit(1)


@cli.command()
@click.option('--data-dir', '-d', required=True,
              type=click.Path(exists=True),
              help='Directory containing training data')
@click.option('--model-name', '-n', required=True,
              help='Name for the trained model')
@click.option('--config', '-c', default='config/setup_config.yaml',
              help='Path to configuration file')
def train(data_dir, model_name, config):
    """Train a new model with custom data."""
    click.echo(f"Training new model '{model_name}' using data from {data_dir}...")
    
    try:
        # Import the training module
        from phytovenomics.training import train_model
        
        # Run training
        result = train_model(data_dir=data_dir, model_name=model_name)
        
        if result:
            click.echo(click.style(f"✅ Model '{model_name}' trained successfully!", fg='green'))
        else:
            click.echo(click.style("❌ Training failed.", fg='red'))
            sys.exit(1)
    except ImportError:
        click.echo(click.style("❌ Could not import training module.", fg='red'))
        click.echo("Make sure Phytovenomics is properly installed.")
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"❌ Training failed with error: {str(e)}", fg='red'))
        sys.exit(1)


def main():
    """Main entry point for the CLI."""
    cli()


if __name__ == "__main__":
    main()