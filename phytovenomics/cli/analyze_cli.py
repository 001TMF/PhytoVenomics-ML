# phytovenomics/cli/analyze_cli.py
#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Phytovenomics Analysis CLI

This module provides a command-line interface for analyzing data with Phytovenomics,
including sequence analysis, compound identification, and visualization tools.
"""

import sys
import click
from pathlib import Path

from phytovenomics.setup.config_manager import ConfigManager


@click.group()
def cli():
    """Phytovenomics analysis commands."""
    pass


@cli.command()
@click.argument('sequence_file', type=click.Path(exists=True))
@click.option('--output', '-o', default='results/analysis',
              help='Directory to save analysis results')
@click.option('--config', '-c', default='config/setup_config.yaml',
              help='Path to configuration file')
def sequence(sequence_file, output, config):
    """Analyze protein or DNA sequence."""
    click.echo(f"Analyzing sequence from {sequence_file}...")
    
    try:
        # Import the sequence analysis module
        from phytovenomics.analysis import analyze_sequence
        
        # Ensure output directory exists
        output_dir = Path(output)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run analysis
        result = analyze_sequence(sequence_file, output_dir=output)
        
        if result:
            click.echo(click.style(f"✅ Sequence analysis completed! Results saved to {output}", fg='green'))
        else:
            click.echo(click.style("❌ Sequence analysis failed.", fg='red'))
            sys.exit(1)
    except ImportError:
        click.echo(click.style("❌ Could not import analysis module.", fg='red'))
        click.echo("Make sure Phytovenomics is properly installed.")
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"❌ Analysis failed with error: {str(e)}", fg='red'))
        sys.exit(1)


@cli.command()
@click.argument('compound_file', type=click.Path(exists=True))
@click.option('--output', '-o', default='results/compounds',
              help='Directory to save compound analysis results')
@click.option('--config', '-c', default='config/setup_config.yaml',
              help='Path to configuration file')
def compounds(compound_file, output, config):
    """Analyze chemical compounds."""
    click.echo(f"Analyzing compounds from {compound_file}...")
    
    try:
        # Import the compound analysis module
        from phytovenomics.analysis import analyze_compounds
        
        # Ensure output directory exists
        output_dir = Path(output)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run analysis
        result = analyze_compounds(compound_file, output_dir=output)
        
        if result:
            click.echo(click.style(f"✅ Compound analysis completed! Results saved to {output}", fg='green'))
        else:
            click.echo(click.style("❌ Compound analysis failed.", fg='red'))
            sys.exit(1)
    except ImportError:
        click.echo(click.style("❌ Could not import analysis module.", fg='red'))
        click.echo("Make sure Phytovenomics is properly installed.")
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"❌ Analysis failed with error: {str(e)}", fg='red'))
        sys.exit(1)


@cli.command()
@click.argument('plant_data', type=click.Path(exists=True))
@click.argument('venom_data', type=click.Path(exists=True))
@click.option('--output', '-o', default='results/interaction',
              help='Directory to save interaction analysis results')
@click.option('--config', '-c', default='config/setup_config.yaml',
              help='Path to configuration file')
def interaction(plant_data, venom_data, output, config):
    """Analyze plant-venom interactions."""
    click.echo(f"Analyzing plant-venom interactions between {plant_data} and {venom_data}...")
    
    try:
        # Import the interaction analysis module
        from phytovenomics.analysis import analyze_interaction
        
        # Ensure output directory exists
        output_dir = Path(output)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run analysis
        result = analyze_interaction(plant_data, venom_data, output_dir=output)
        
        if result:
            click.echo(click.style(f"✅ Interaction analysis completed! Results saved to {output}", fg='green'))
        else:
            click.echo(click.style("❌ Interaction analysis failed.", fg='red'))
            sys.exit(1)
    except ImportError:
        click.echo(click.style("❌ Could not import analysis module.", fg='red'))
        click.echo("Make sure Phytovenomics is properly installed.")
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"❌ Analysis failed with error: {str(e)}", fg='red'))
        sys.exit(1)


@cli.command()
@click.argument('data_file', type=click.Path(exists=True))
@click.option('--type', '-t', required=True, 
              type=click.Choice(['sequence', 'compound', 'interaction']),
              help='Type of data to visualize')
@click.option('--output', '-o', default='results/visualizations',
              help='Directory to save visualizations')
@click.option('--format', '-f', default='png',
              type=click.Choice(['png', 'svg', 'pdf']),
              help='Output format for visualizations')
@click.option('--config', '-c', default='config/setup_config.yaml',
              help='Path to configuration file')
def visualize(data_file, type, output, format, config):
    """Generate visualizations from analysis data."""
    click.echo(f"Generating {type} visualizations from {data_file}...")
    
    try:
        # Import the visualization module
        from phytovenomics.visualization import create_visualization
        
        # Ensure output directory exists
        output_dir = Path(output)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run visualization
        result = create_visualization(data_file, vis_type=type, 
                                     output_dir=output, output_format=format)
        
        if result:
            click.echo(click.style(f"✅ Visualizations created! Results saved to {output}", fg='green'))
        else:
            click.echo(click.style("❌ Visualization creation failed.", fg='red'))
            sys.exit(1)
    except ImportError:
        click.echo(click.style("❌ Could not import visualization module.", fg='red'))
        click.echo("Make sure Phytovenomics is properly installed.")
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"❌ Visualization failed with error: {str(e)}", fg='red'))
        sys.exit(1)


def main():
    """Main entry point for the CLI."""
    cli()


if __name__ == "__main__":
    main()