#!/usr/bin/env python
# utils/validation_metrics.py

"""
Validation metrics tracker for the Phytovenomics platform.
This module provides utilities for tracking, analyzing, and visualizing metrics
from the various pipeline components.
"""

import os
import json
import logging
import time
from typing import Dict, List, Any, Optional, Union, Tuple
import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)

class ValidationMetrics:
    """
    Class for tracking validation metrics across the pipeline.
    Allows adding, retrieving, and visualizing various metrics.
    """
    
    def __init__(self, output_dir: Optional[str] = None):
        """
        Initialize the metrics tracker.
        
        Args:
            output_dir: Directory to save metrics and visualizations
        """
        self.metrics = {}
        self.output_dir = output_dir
        
        # Create output directory if provided
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        
        # Initialize timestamp
        self.start_time = time.time()
        self.timestamps = {}
    
    def add_metric(self, name: str, value: Any) -> None:
        """
        Add a metric value to tracking.
        
        Args:
            name: Name of the metric
            value: Value of the metric (number, string, boolean, etc.)
        """
        # If metric doesn't exist yet, initialize it as a list
        if name not in self.metrics:
            self.metrics[name] = []
            self.timestamps[name] = []
        
        # Add the value and timestamp
        self.metrics[name].append(value)
        self.timestamps[name].append(time.time() - self.start_time)
        
        logger.debug(f"Added metric: {name} = {value}")
    
    def get_metric(self, name: str) -> List[Any]:
        """
        Get all values for a metric.
        
        Args:
            name: Name of the metric
            
        Returns:
            List of all values for the metric
        """
        return self.metrics.get(name, [])
    
    def get_metric_average(self, name: str) -> float:
        """
        Get the average value of a numeric metric.
        
        Args:
            name: Name of the metric
            
        Returns:
            Average value or 0 if no values or non-numeric
        """
        values = self.get_metric(name)
        
        if not values:
            return 0
        
        try:
            return sum(values) / len(values)
        except (TypeError, ValueError):
            logger.warning(f"Cannot calculate average for non-numeric metric: {name}")
            return 0
    
    def get_metric_latest(self, name: str) -> Any:
        """
        Get the most recent value for a metric.
        
        Args:
            name: Name of the metric
            
        Returns:
            Latest value or None if metric doesn't exist
        """
        values = self.get_metric(name)
        
        if not values:
            return None
        
        return values[-1]
    
    def plot_metrics(self, metric_names: List[str], output_path: Optional[str] = None,
                    title: str = "Metrics", xlabel: str = "Metric", 
                    ylabel: str = "Value") -> str:
        """
        Create a bar plot for multiple metrics.
        
        Args:
            metric_names: List of metric names to plot
            output_path: Path to save the plot
            title: Plot title
            xlabel: X-axis label
            ylabel: Y-axis label
            
        Returns:
            Path to the saved plot or empty string if not saved
        """
        # Collect values for each metric (use average if multiple values)
        values = []
        labels = []
        
        for name in metric_names:
            if name in self.metrics:
                values.append(self.get_metric_average(name))
                # Format the label for display
                display_name = name.replace('_', ' ').title()
                labels.append(display_name)
        
        # Create the plot if we have data
        if not values:
            logger.warning("No metrics to plot")
            return ""
        
        plt.figure(figsize=(10, 6))
        bars = plt.bar(labels, values)
        
        # Add value labels on top of bars
        for bar in bars:
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width() / 2,
                height + 0.01 * max(values),
                f'{height:.2f}',
                ha='center',
                fontsize=9
            )
        
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        # Save plot if output path provided
        if output_path:
            # Create directory if needed
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            plt.savefig(output_path)
            logger.info(f"Saved metrics plot to {output_path}")
            plt.close()
            return output_path
        
        return ""
    
    def plot_time_series(self, metric_name: str, output_path: Optional[str] = None,
                        title: Optional[str] = None, xlabel: str = "Time (s)", 
                        ylabel: Optional[str] = None) -> str:
        """
        Create a time series plot for a metric.
        
        Args:
            metric_name: Name of the metric to plot
            output_path: Path to save the plot
            title: Plot title (defaults to metric name)
            xlabel: X-axis label
            ylabel: Y-axis label (defaults to metric name)
            
        Returns:
            Path to the saved plot or empty string if not saved
        """
        if metric_name not in self.metrics:
            logger.warning(f"Metric {metric_name} not found for time series plot")
            return ""
        
        values = self.metrics[metric_name]
        times = self.timestamps[metric_name]
        
        # Use defaults if not provided
        if title is None:
            title = f"{metric_name.replace('_', ' ').title()} Over Time"
        
        if ylabel is None:
            ylabel = metric_name.replace('_', ' ').title()
        
        plt.figure(figsize=(10, 6))
        plt.plot(times, values, 'o-')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        # Save plot if output path provided
        if output_path:
            # Create directory if needed
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            plt.savefig(output_path)
            logger.info(f"Saved time series plot to {output_path}")
            plt.close()
            return output_path
        
        return ""
    
    def export_to_json(self, output_path: str) -> Dict[str, Any]:
        """
        Export all metrics to a JSON file.
        
        Args:
            output_path: Path to save the JSON file
            
        Returns:
            Dictionary of metrics
        """
        # Prepare the metrics dictionary with summary statistics
        metrics_export = {}
        
        for name, values in self.metrics.items():
            # Skip empty metrics
            if not values:
                continue
            
            # Try to calculate statistics for numeric metrics
            try:
                avg = sum(values) / len(values)
                min_val = min(values)
                max_val = max(values)
                
                metrics_export[name] = {
                    "values": values,
                    "timestamps": self.timestamps[name],
                    "summary": {
                        "average": avg,
                        "min": min_val,
                        "max": max_val,
                        "count": len(values),
                        "latest": values[-1]
                    }
                }
            except (TypeError, ValueError):
                # For non-numeric metrics, just store the values
                metrics_export[name] = {
                    "values": values,
                    "timestamps": self.timestamps[name],
                    "summary": {
                        "count": len(values),
                        "latest": values[-1]
                    }
                }
        
        # Create directory if needed
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Save to file
        with open(output_path, 'w') as f:
            json.dump(metrics_export, f, indent=2)
        
        logger.info(f"Exported metrics to {output_path}")
        return metrics_export
    
    def summary(self) -> Dict[str, Any]:
        """
        Generate a summary of all metrics.
        
        Returns:
            Dictionary with metric summaries
        """
        summary_dict = {}
        
        for name, values in self.metrics.items():
            # Skip empty metrics
            if not values:
                continue
            
            # Try to calculate statistics for numeric metrics
            try:
                avg = sum(values) / len(values)
                min_val = min(values)
                max_val = max(values)
                
                summary_dict[name] = {
                    "average": avg,
                    "min": min_val,
                    "max": max_val,
                    "count": len(values),
                    "latest": values[-1]
                }
            except (TypeError, ValueError):
                # For non-numeric metrics, just store count and latest
                summary_dict[name] = {
                    "count": len(values),
                    "latest": values[-1]
                }
        
        return summary_dict
    
    def print_summary(self) -> None:
        """
        Print a summary of all metrics to the console.
        """
        summary = self.summary()
        
        print("\n=== Metrics Summary ===")
        
        for name, stats in summary.items():
            print(f"\n{name.replace('_', ' ').title()}:")
            
            for stat_name, value in stats.items():
                if isinstance(value, (int, float)) and stat_name != "count":
                    print(f"  {stat_name.title()}: {value:.4f}")
                else:
                    print(f"  {stat_name.title()}: {value}")
        
        print("\n=====================")


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)
    
    # Create metrics tracker
    metrics = ValidationMetrics("test_metrics")
    
    # Add some test metrics
    for i in range(10):
        metrics.add_metric("accuracy", 0.75 + i * 0.02)
        metrics.add_metric("loss", 0.5 - i * 0.03)
        metrics.add_metric("iteration", i + 1)
        time.sleep(0.1)  # Simulate time passing
    
    # Plot metrics
    metrics.plot_metrics(["accuracy", "loss"], "test_metrics/metrics_plot.png")
    metrics.plot_time_series("accuracy", "test_metrics/accuracy_time_series.png")
    
    # Export metrics
    metrics.export_to_json("test_metrics/metrics.json")
    
    # Print summary
    metrics.print_summary()