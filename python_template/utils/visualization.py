#!/usr/bin/env python
# utils/visualization.py

import os
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional, Union, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)

class PlotHelper:
    """
    Helper class for creating common plots for visualization.
    """
    
    @staticmethod
    def plot_bar_chart(data, labels, title="", xlabel="", ylabel="", output_path=None, figsize=(10, 6), color=None):
        """
        Create a bar chart.
        
        Args:
            data: Values for the bars
            labels: Labels for the bars
            title: Title of the plot
            xlabel: Label for the x-axis
            ylabel: Label for the y-axis
            output_path: Path to save the plot
            figsize: Figure size (width, height)
            color: Bar color
        """
        plt.figure(figsize=figsize)
        bars = plt.bar(labels, data, color=color)
        
        # Add labels and title
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        
        # Add value labels on top of bars
        for bar in bars:
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width() / 2,
                height + 0.01 * max(data),
                f'{height:.2f}',
                ha='center',
                fontsize=9
            )
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path)
            logger.info(f"Saved bar chart to {output_path}")
            plt.close()
            return output_path
        
        return plt.gcf()
    
    @staticmethod
    def plot_line_chart(data, x_values=None, labels=None, title="", xlabel="", ylabel="", 
                       output_path=None, figsize=(10, 6), marker='o', linewidth=2):
        """
        Create a line chart.
        
        Args:
            data: List of data series
            x_values: X-axis values
            labels: Labels for each data series
            title: Title of the plot
            xlabel: Label for the x-axis
            ylabel: Label for the y-axis
            output_path: Path to save the plot
            figsize: Figure size (width, height)
            marker: Marker style
            linewidth: Width of the line
        """
        plt.figure(figsize=figsize)
        
        # Handle single data series
        if not isinstance(data[0], (list, tuple, np.ndarray)):
            data = [data]
            labels = labels or [""]
        
        # Create x values if not provided
        if x_values is None:
            x_values = list(range(len(data[0])))
        
        # Plot each data series
        for i, y_values in enumerate(data):
            label = labels[i] if labels and i < len(labels) else f"Series {i+1}"
            plt.plot(x_values, y_values, marker=marker, linewidth=linewidth, label=label)
        
        if labels and len(labels) > 1:
            plt.legend()
        
        # Add labels and title
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True, linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path)
            logger.info(f"Saved line chart to {output_path}")
            plt.close()
            return output_path
        
        return plt.gcf()
    
    @staticmethod
    def plot_scatter(x_data, y_data, labels=None, colors=None, sizes=None, 
                    title="", xlabel="", ylabel="", output_path=None, figsize=(10, 6)):
        """
        Create a scatter plot.
        
        Args:
            x_data: X coordinates
            y_data: Y coordinates
            labels: Labels for the scatter points
            colors: Colors for the scatter points
            sizes: Sizes for the scatter points
            title: Title of the plot
            xlabel: Label for the x-axis
            ylabel: Label for the y-axis
            output_path: Path to save the plot
            figsize: Figure size (width, height)
        """
        plt.figure(figsize=figsize)
        
        scatter = plt.scatter(
            x_data, y_data, 
            c=colors, 
            s=sizes or 50, 
            alpha=0.7, 
            edgecolors='w'
        )
        
        # Add labels if provided
        if labels is not None:
            for i, label in enumerate(labels):
                plt.annotate(
                    label,
                    (x_data[i], y_data[i]),
                    textcoords="offset points",
                    xytext=(0, 5),
                    ha='center'
                )
        
        # Add labels and title
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True, linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path)
            logger.info(f"Saved scatter plot to {output_path}")
            plt.close()
            return output_path
        
        return plt.gcf()
    
    @staticmethod
    def plot_heatmap(data, x_labels=None, y_labels=None, title="", xlabel="", ylabel="", 
                    output_path=None, figsize=(10, 8), cmap="viridis", annot=True):
        """
        Create a heatmap.
        
        Args:
            data: 2D array of values
            x_labels: Labels for the x-axis
            y_labels: Labels for the y-axis
            title: Title of the plot
            xlabel: Label for the x-axis
            ylabel: Label for the y-axis
            output_path: Path to save the plot
            figsize: Figure size (width, height)
            cmap: Colormap for the heatmap
            annot: Whether to annotate cells
        """
        plt.figure(figsize=figsize)
        
        # Create heatmap
        plt.imshow(data, cmap=cmap)
        
        # Add colorbar
        plt.colorbar(label='Value')
        
        # Set x and y labels if provided
        if x_labels is not None:
            plt.xticks(range(len(x_labels)), x_labels, rotation=45, ha="right")
        if y_labels is not None:
            plt.yticks(range(len(y_labels)), y_labels)
        
        # Add annotations if requested
        if annot:
            for i in range(len(data)):
                for j in range(len(data[i])):
                    plt.text(j, i, f"{data[i][j]:.2f}", ha="center", va="center", color="w")
        
        # Add labels and title
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path)
            logger.info(f"Saved heatmap to {output_path}")
            plt.close()
            return output_path
        
        return plt.gcf()

    @staticmethod
    def plot_3d_scatter(x_data, y_data, z_data, labels=None, colors=None, sizes=None, 
                      title="", xlabel="", ylabel="", zlabel="", output_path=None, figsize=(10, 8)):
        """
        Create a 3D scatter plot.
        
        Args:
            x_data: X coordinates
            y_data: Y coordinates
            z_data: Z coordinates
            labels: Labels for the scatter points
            colors: Colors for the scatter points
            sizes: Sizes for the scatter points
            title: Title of the plot
            xlabel: Label for the x-axis
            ylabel: Label for the y-axis
            zlabel: Label for the z-axis
            output_path: Path to save the plot
            figsize: Figure size (width, height)
        """
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        
        # Create scatter plot
        scatter = ax.scatter(
            x_data, y_data, z_data,
            c=colors,
            s=sizes or 50,
            alpha=0.7
        )
        
        # Add labels if provided
        if labels is not None:
            for i, label in enumerate(labels):
                ax.text(x_data[i], y_data[i], z_data[i], label)
        
        # Add labels and title
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        ax.set_title(title)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path)
            logger.info(f"Saved 3D scatter plot to {output_path}")
            plt.close()
            return output_path
        
        return fig


class AntibodyVisualizer:
    """
    Class for visualizing antibody structures and properties.
    """
    
    def __init__(self, output_dir: Optional[str] = None):
        """
        Initialize the antibody visualizer.
        
        Args:
            output_dir: Directory to save visualizations
        """
        self.output_dir = output_dir
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
    
    def plot_binding_energy(self, antibodies, target_name=None, output_path=None):
        """
        Plot binding energy of antibodies against a target.
        
        Args:
            antibodies: List of antibody objects or dictionary with binding energies
            target_name: Name of the target antigen/epitope
            output_path: Path to save the plot
        """
        # Extract binding energies and names
        if hasattr(antibodies[0], 'binding_energy'):
            names = [ab.name for ab in antibodies]
            energies = [ab.binding_energy for ab in antibodies]
        else:
            names = [ab.get('name', f"Ab{i}") for i, ab in enumerate(antibodies)]
            energies = [ab.get('binding_energy', 0) for ab in antibodies]
        
        # Sort by binding energy (stronger binding first)
        sorted_indices = np.argsort(energies)
        sorted_names = [names[i] for i in sorted_indices]
        sorted_energies = [energies[i] for i in sorted_indices]
        
        # Create the plot
        if output_path is None and self.output_dir:
            output_path = os.path.join(self.output_dir, 'binding_energy.png')
            
        target_str = f" against {target_name}" if target_name else ""
        return PlotHelper.plot_bar_chart(
            sorted_energies, 
            sorted_names,
            title=f"Binding Energy{target_str}",
            xlabel="Antibody",
            ylabel="Energy (kcal/mol)",
            output_path=output_path,
            figsize=(12, 6),
            color='skyblue'
        )

    def plot_antibody_metrics(self, antibodies, metrics=None, output_path=None):
        """
        Plot multiple metrics for a set of antibodies.
        
        Args:
            antibodies: List of antibody objects
            metrics: List of metrics to plot
            output_path: Path to save the plot
        """
        if metrics is None:
            metrics = ['binding_energy', 'stability', 'developability', 'manufacturability']
        
        # Prepare data
        ab_names = []
        data_by_metric = {metric: [] for metric in metrics}
        
        for ab in antibodies:
            if hasattr(ab, 'name'):
                ab_names.append(ab.name)
                for metric in metrics:
                    data_by_metric[metric].append(getattr(ab, metric, 0))
            else:
                ab_names.append(ab.get('name', f"Ab{len(ab_names)+1}"))
                for metric in metrics:
                    data_by_metric[metric].append(ab.get(metric, 0))
        
        # Set up plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Set width of bars
        bar_width = 0.2
        positions = np.arange(len(ab_names))
        
        # Plot bars for each metric
        for i, (metric, values) in enumerate(data_by_metric.items()):
            offset = (i - len(metrics)/2 + 0.5) * bar_width
            ax.bar(positions + offset, values, bar_width, label=metric.replace('_', ' ').title())
        
        # Add labels, title and legend
        ax.set_xlabel('Antibody')
        ax.set_ylabel('Value')
        ax.set_title('Antibody Properties Comparison')
        ax.set_xticks(positions)
        ax.set_xticklabels(ab_names, rotation=45, ha='right')
        ax.legend()
        
        plt.tight_layout()
        
        if output_path is None and self.output_dir:
            output_path = os.path.join(self.output_dir, 'antibody_metrics.png')
            
        if output_path:
            plt.savefig(output_path)
            logger.info(f"Saved antibody metrics plot to {output_path}")
            plt.close()
            return output_path
        
        return fig

    def plot_evolution_progress(self, generations, metric='fitness', output_path=None):
        """
        Plot the progress of evolutionary optimization.
        
        Args:
            generations: List of generation data
            metric: Metric to plot
            output_path: Path to save the plot
        """
        gen_numbers = []
        best_values = []
        avg_values = []
        
        for i, gen in enumerate(generations):
            gen_numbers.append(i + 1)
            
            if isinstance(gen, list):
                # Extract values from list of individuals
                values = [getattr(ind, metric, 0) if hasattr(ind, metric) else ind.get(metric, 0) 
                         for ind in gen]
                best_values.append(max(values))
                avg_values.append(sum(values) / len(values) if values else 0)
            else:
                # Data is already processed
                best_values.append(gen.get(f'best_{metric}', 0))
                avg_values.append(gen.get(f'avg_{metric}', 0))
        
        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.plot(gen_numbers, best_values, 'o-', label=f'Best {metric}', linewidth=2)
        plt.plot(gen_numbers, avg_values, 'o--', label=f'Average {metric}', linewidth=2)
        
        plt.xlabel('Generation')
        plt.ylabel(metric.replace('_', ' ').title())
        plt.title(f'Evolution of {metric.replace("_", " ").title()} Over Generations')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        
        if output_path is None and self.output_dir:
            output_path = os.path.join(self.output_dir, f'evolution_{metric}.png')
            
        if output_path:
            plt.savefig(output_path)
            logger.info(f"Saved evolution progress plot to {output_path}")
            plt.close()
            return output_path
        
        return plt.gcf()
    
    def plot_cdr_profile(self, antibody, output_path=None):
        """
        Plot CDR properties for an antibody.
        
        Args:
            antibody: Antibody object or dictionary
            output_path: Path to save the plot
        """
        # Calculate properties for CDRs
        cdr_names = ["H1", "H2", "H3", "L1", "L2", "L3"]
        
        # Mock CDR property calculation
        if hasattr(antibody, 'cdr_sequences'):
            cdr_sequences = antibody.cdr_sequences
        else:
            cdr_sequences = antibody.get('cdr_sequences', {})
            
        # Calculate properties
        hydrophobicity = []
        charge = []
        length = []
        
        # Mock calculation of properties
        for cdr in cdr_names:
            if cdr in cdr_sequences:
                seq = cdr_sequences[cdr]
                length.append(len(seq))
                # Mock calculations
                hydrophobicity.append(0.5 + 0.3 * (hash(seq) % 10) / 10)
                charge.append(-2 + (hash(seq) % 5))
            else:
                length.append(0)
                hydrophobicity.append(0)
                charge.append(0)
        
        # Create the plot
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
        
        # Plot hydrophobicity
        ax1.bar(cdr_names, hydrophobicity, color='skyblue')
        ax1.set_ylabel('Hydrophobicity')
        ax1.set_title('CDR Properties')
        
        # Plot charge
        ax2.bar(cdr_names, charge, color='salmon')
        ax2.set_ylabel('Net Charge')
        
        # Plot length
        ax3.bar(cdr_names, length, color='lightgreen')
        ax3.set_ylabel('Length')
        ax3.set_xlabel('CDR')
        
        plt.tight_layout()
        
        if output_path is None and self.output_dir:
            if hasattr(antibody, 'name'):
                name = antibody.name
            else:
                name = antibody.get('name', 'antibody')
            output_path = os.path.join(self.output_dir, f'{name}_cdr_profile.png')
            
        if output_path:
            plt.savefig(output_path)
            logger.info(f"Saved CDR profile plot to {output_path}")
            plt.close()
            return output_path
        
        return fig


class DataVisualizer:
    """
    Class for visualizing general data analysis and metrics.
    """
    
    def __init__(self, output_dir: Optional[str] = None):
        """
        Initialize the data visualizer.
        
        Args:
            output_dir: Directory to save visualizations
        """
        self.output_dir = output_dir
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
    
    def plot_metrics_comparison(self, metrics_dict, title="Metrics Comparison", 
                              output_path=None, figsize=(12, 8)):
        """
        Create a comparison plot of multiple metrics.
        
        Args:
            metrics_dict: Dictionary of metrics {name: value}
            title: Plot title
            output_path: Path to save the plot
            figsize: Figure size (width, height)
        """
        names = list(metrics_dict.keys())
        values = list(metrics_dict.values())
        
        # Create the plot
        plt.figure(figsize=figsize)
        bars = plt.barh(names, values, color='skyblue')
        
        # Add value annotations
        for bar in bars:
            width = bar.get_width()
            plt.text(width + 0.01 * max(values), 
                    bar.get_y() + bar.get_height()/2, 
                    f'{width:.2f}', 
                    va='center')
        
        plt.xlabel('Value')
        plt.title(title)
        plt.grid(True, linestyle='--', alpha=0.7, axis='x')
        
        plt.tight_layout()
        
        if output_path is None and self.output_dir:
            safe_title = title.replace(' ', '_').lower()
            output_path = os.path.join(self.output_dir, f'{safe_title}.png')
            
        if output_path:
            plt.savefig(output_path)
            logger.info(f"Saved metrics comparison to {output_path}")
            plt.close()
            return output_path
        
        return plt.gcf()
    
    def plot_time_series(self, time_points, data_series, labels=None, title="Time Series", 
                       xlabel="Time", ylabel="Value", output_path=None, figsize=(12, 6)):
        """
        Create a time series plot.
        
        Args:
            time_points: List of time points
            data_series: List of data series (list of lists)
            labels: Labels for each data series
            title: Plot title
            xlabel: X-axis label
            ylabel: Y-axis label
            output_path: Path to save the plot
            figsize: Figure size (width, height)
        """
        return PlotHelper.plot_line_chart(
            data=data_series,
            x_values=time_points,
            labels=labels,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
            output_path=output_path,
            figsize=figsize
        )
    
    def create_report(self, metrics, output_path=None, title="Analysis Report"):
        """
        Create a simple markdown report from metrics.
        
        Args:
            metrics: Dictionary of metrics
            output_path: Path to save the report
            title: Report title
        """
        report = [f"# {title}", ""]
        
        # Add timestamp
        from datetime import datetime
        report.append(f"Report generated at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report.append("")
        
        # Add sections for different metric categories
        categories = {}
        for key, value in metrics.items():
            # Group by prefix (before first underscore)
            if '_' in key:
                category = key.split('_')[0]
            else:
                category = 'general'
                
            if category not in categories:
                categories[category] = []
                
            categories[category].append((key, value))
        
        # Write each category
        for category, items in categories.items():
            report.append(f"## {category.title()}")
            report.append("")
            
            for key, value in items:
                display_name = key.replace(f"{category}_", "", 1).replace('_', ' ').title()
                if isinstance(value, (int, float)):
                    report.append(f"- {display_name}: {value:.4f}")
                else:
                    report.append(f"- {display_name}: {value}")
            
            report.append("")
        
        # Write to file if output path provided
        if output_path is None and self.output_dir:
            output_path = os.path.join(self.output_dir, 'report.md')
            
        if output_path:
            with open(output_path, 'w') as f:
                f.write('\n'.join(report))
            logger.info(f"Saved report to {output_path}")
            return output_path
        
        return '\n'.join(report)


def visualize_antibody_binding(antibody_list, toxin_list, output_path=None):
    """
    Visualize binding between antibodies and toxins.
    
    Args:
        antibody_list: List of antibodies
        toxin_list: List of toxins
        output_path: Path to save visualization
    
    Returns:
        Path to saved visualization
    """
    # Create binding matrix
    ab_names = [ab['id'] if isinstance(ab, dict) else ab.name for ab in antibody_list]
    toxin_names = [t['id'] if isinstance(t, dict) else t.name for t in toxin_list]
    
    binding_matrix = np.zeros((len(ab_names), len(toxin_names)))
    
    # Fill binding matrix with mock binding strengths
    for i, ab in enumerate(antibody_list):
        for j, toxin in enumerate(toxin_list):
            # Generate consistent mock binding strength
            id_str = f"{ab_names[i]}_{toxin_names[j]}"
            binding_matrix[i, j] = 0.3 + 0.6 * (hash(id_str) % 100) / 100
    
    # Create heatmap
    plt.figure(figsize=(12, 8))
    plt.imshow(binding_matrix, cmap='viridis')
    plt.colorbar(label='Binding Strength')
    
    # Set labels
    plt.xticks(range(len(toxin_names)), toxin_names, rotation=45, ha='right')
    plt.yticks(range(len(ab_names)), ab_names)
    
    plt.xlabel('Toxins')
    plt.ylabel('Antibodies')
    plt.title('Antibody-Toxin Binding Matrix')
    
    plt.tight_layout()
    
    # Save if output path is provided
    if output_path:
        plt.savefig(output_path)
        logger.info(f"Saved binding matrix visualization to {output_path}")
        plt.close()
        return output_path
    
    return plt.gcf()


def visualize_antibody_structure(pdb_data=None, highlight_regions=None, output_path=None):
    """
    Create a mock visualization of antibody structure.
    
    Args:
        pdb_data: PDB data or file path
        highlight_regions: Dictionary mapping region names to residue ranges
        output_path: Path to save visualization
    
    Returns:
        Path to saved visualization or matplotlib figure
    """
    # Create a mock 3D structure plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Generate random coordinates for a mock antibody structure
    np.random.seed(42)  # For reproducibility
    n_points = 100
    x = np.random.normal(0, 1, n_points)
    y = np.random.normal(0, 1, n_points)
    z = np.random.normal(0, 1, n_points)
    
    # Plot the main structure
    ax.plot(x, y, z, 'o-', color='gray', alpha=0.6, linewidth=1, markersize=3)
    
    # Highlight specific regions if provided
    if highlight_regions:
        colors = ['red', 'green', 'blue', 'orange', 'purple', 'cyan']
        for i, (region_name, residue_range) in enumerate(highlight_regions.items()):
            color = colors[i % len(colors)]
            start, end = residue_range
            
            # Generate coordinates for highlighted region
            region_x = x[start:end] + 0.1
            region_y = y[start:end] + 0.1
            region_z = z[start:end] + 0.1
            
            ax.plot(region_x, region_y, region_z, 'o-', 
                  color=color, 
                  linewidth=2, 
                  markersize=5, 
                  label=region_name)
    
    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Antibody Structure')
    
    # Add legend if regions are highlighted
    if highlight_regions:
        ax.legend()
    
    # Set the viewing angle
    ax.view_init(elev=30, azim=45)
    
    # Save if output path is provided
    if output_path:
        plt.savefig(output_path)
        logger.info(f"Saved structure visualization to {output_path}")
        plt.close()
        return output_path
    
    return fig


# For backwards compatibility
def plot_metrics(metrics_dict, title="Metrics", output_path=None):
    """
    Plot metrics from a dictionary.
    
    Args:
        metrics_dict: Dictionary of metrics
        title: Plot title
        output_path: Path to save the plot
    
    Returns:
        Path to saved plot or matplotlib figure
    """
    visualizer = DataVisualizer()
    return visualizer.plot_metrics_comparison(metrics_dict, title, output_path)