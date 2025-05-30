#!/usr/bin/env python
# utils/validation_utils.py

import os
import json
import time
import logging
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Union, Any, Tuple
from pathlib import Path
from datetime import datetime
import pandas as pd
from collections import defaultdict
from scipy import stats

logger = logging.getLogger("Phytovenomics.ValidationUtils")

class ValidationMetricsTracker:
    """
    Comprehensive validation metrics tracking for the ML platform.
    
    This class enables tracking, recording, and visualizing various metrics:
    1. Binding affinity predictions
    2. Structure prediction quality
    3. Antibody design quality
    4. Toxin neutralization efficacy
    5. Model performance over time
    """
    
    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize the metrics tracker with configuration.
        
        Args:
            config: Optional configuration dictionary
        """
        self.config = config or {}
        
        # Default configurations
        self.metrics_dir = self.config.get("metrics_dir", "output/metrics")
        self.metrics_history = defaultdict(list)
        self.current_run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Ensure output directory exists
        os.makedirs(self.metrics_dir, exist_ok=True)
        
        logger.info(f"Initialized validation metrics tracker with run ID: {self.current_run_id}")
    
    def calculate_binding_metrics(self, predicted_affinities: List[float], 
                                 true_affinities: List[float]) -> Dict[str, float]:
        """
        Calculate metrics for binding affinity prediction quality.
        
        Args:
            predicted_affinities: List of predicted binding affinities (in nM)
            true_affinities: List of true binding affinities (in nM)
            
        Returns:
            Dictionary of metrics
        """
        if len(predicted_affinities) != len(true_affinities) or len(predicted_affinities) == 0:
            logger.error("Predicted and true affinity lists must be non-empty and of the same length")
            return {}
            
        # Convert to log space for better comparison (binding affinities are typically log-distributed)
        log_pred = np.log10(np.array(predicted_affinities) + 1e-10)
        log_true = np.log10(np.array(true_affinities) + 1e-10)
        
        # Calculate metrics
        metrics = {
            "mae": float(np.mean(np.abs(np.array(predicted_affinities) - np.array(true_affinities)))),
            "rmse": float(np.sqrt(np.mean((np.array(predicted_affinities) - np.array(true_affinities))**2))),
            "pearson_r": float(stats.pearsonr(predicted_affinities, true_affinities)[0]),
            "spearman_rho": float(stats.spearmanr(predicted_affinities, true_affinities)[0]),
            "log_mae": float(np.mean(np.abs(log_pred - log_true))),
            "log_rmse": float(np.sqrt(np.mean((log_pred - log_true)**2))),
            "hits_at_10": self._calculate_hits_at_k(predicted_affinities, true_affinities, 10),
        }
        
        # Add to metrics history
        self._add_to_history("binding_metrics", metrics)
        
        return metrics
    
    def calculate_structure_metrics(self, predicted_structures: List[Dict], 
                                  reference_structures: List[Dict]) -> Dict[str, float]:
        """
        Calculate metrics for structure prediction quality.
        
        Args:
            predicted_structures: List of predicted structure dictionaries
            reference_structures: List of reference structure dictionaries
            
        Returns:
            Dictionary of metrics
        """
        if len(predicted_structures) != len(reference_structures) or len(predicted_structures) == 0:
            logger.error("Predicted and reference structure lists must be non-empty and of the same length")
            return {}
            
        # Extract metrics from structures (assuming each has required fields)
        # This is a simplified implementation - actual code would compute RMSD, TM-score, etc.
        try:
            # Extract confidence scores
            confidence_scores = [s.get("confidence", 0) for s in predicted_structures]
            
            # Extract pLDDT scores if available
            plddt_scores = []
            for structure in predicted_structures:
                if "plddt" in structure:
                    if isinstance(structure["plddt"], list):
                        plddt_scores.extend(structure["plddt"])
                    else:
                        plddt_scores.append(structure["plddt"])
                        
            metrics = {
                "avg_confidence": float(np.mean(confidence_scores)),
                "min_confidence": float(np.min(confidence_scores)),
                "max_confidence": float(np.max(confidence_scores)),
            }
            
            if plddt_scores:
                metrics["avg_plddt"] = float(np.mean(plddt_scores))
                metrics["plddt_below_50"] = float(np.mean(np.array(plddt_scores) < 50))
                
            # In a real implementation, we would calculate RMSD between predicted and reference structures
            # metrics["avg_rmsd"] = calculate_average_rmsd(predicted_structures, reference_structures)
            # metrics["avg_tm_score"] = calculate_average_tm_score(predicted_structures, reference_structures)
            
            # Add to metrics history
            self._add_to_history("structure_metrics", metrics)
            
            return metrics
            
        except Exception as e:
            logger.error(f"Error calculating structure metrics: {str(e)}")
            return {"error": str(e)}
    
    def calculate_antibody_design_metrics(self, 
                                        designed_antibodies: List[Dict], 
                                        target_properties: Optional[Dict] = None) -> Dict[str, float]:
        """
        Calculate metrics for antibody design quality.
        
        Args:
            designed_antibodies: List of designed antibody dictionaries
            target_properties: Optional target property specifications
            
        Returns:
            Dictionary of metrics
        """
        if not designed_antibodies:
            logger.error("No designed antibodies provided")
            return {}
        
        # Extract developability scores
        developability_scores = [ab.get("developability_score", 0) for ab in designed_antibodies]
        
        # Extract binding affinities for all toxins
        binding_affinities = []
        for ab in designed_antibodies:
            binds_to = ab.get("binds_to", [])
            for binding in binds_to:
                if "affinity_nm" in binding:
                    binding_affinities.append(binding["affinity_nm"])
        
        # Calculate metrics
        metrics = {
            "num_designs": len(designed_antibodies),
            "avg_developability": float(np.mean(developability_scores)) if developability_scores else 0,
        }
        
        if binding_affinities:
            metrics["avg_affinity"] = float(np.mean(binding_affinities))
            metrics["min_affinity"] = float(np.min(binding_affinities))  # Lower is better
            metrics["max_affinity"] = float(np.max(binding_affinities))
            
        # Calculate CDR length distributions
        cdr_lengths = defaultdict(list)
        for ab in designed_antibodies:
            if "cdr_lengths" in ab:
                for cdr, length in ab["cdr_lengths"].items():
                    cdr_lengths[cdr].append(length)
        
        for cdr, lengths in cdr_lengths.items():
            metrics[f"{cdr}_avg_length"] = float(np.mean(lengths))
        
        # Calculate sequence novelty (if there are target properties)
        if target_properties and "reference_sequences" in target_properties:
            novelty_scores = []
            for ab in designed_antibodies:
                if "sequence" in ab and ab["sequence"]:
                    max_similarity = self._max_sequence_similarity(
                        ab["sequence"], target_properties["reference_sequences"]
                    )
                    # Novelty is inverse of similarity
                    novelty_scores.append(1.0 - max_similarity)
                    
            if novelty_scores:
                metrics["avg_novelty"] = float(np.mean(novelty_scores))
                metrics["min_novelty"] = float(np.min(novelty_scores))
        
        # Add to metrics history
        self._add_to_history("antibody_design_metrics", metrics)
        
        return metrics
    
    def calculate_neutralization_metrics(self, 
                                       cocktail: Dict, 
                                       toxin_families: Dict[str, List[str]]) -> Dict[str, float]:
        """
        Calculate metrics for toxin neutralization efficacy.
        
        Args:
            cocktail: Dictionary containing antibody cocktail information
            toxin_families: Dictionary mapping toxin family names to toxin IDs
            
        Returns:
            Dictionary of metrics
        """
        if not cocktail or "coverage_metrics" not in cocktail:
            logger.error("Invalid cocktail data provided")
            return {}
        
        # Extract metrics from cocktail
        coverage_metrics = cocktail.get("coverage_metrics", {})
        
        # Calculate overall metrics
        metrics = {
            "antibody_count": cocktail.get("antibody_count", 0),
            "average_coverage": cocktail.get("average_coverage", 0),
            "manufacturability_score": cocktail.get("manufacturability_score", 0),
        }
        
        # Add family-specific coverage
        family_coverage = coverage_metrics.get("family_coverage", {})
        for family, coverage in family_coverage.items():
            metrics[f"family_coverage_{family}"] = coverage
        
        # Calculate cost-effectiveness (simplified)
        if metrics["antibody_count"] > 0:
            # Assume cost scales with number of antibodies
            estimated_cost = metrics["antibody_count"] * 100  # Arbitrary cost units
            metrics["cost_per_coverage_point"] = estimated_cost / (metrics["average_coverage"] * 100)
        
        # Add to metrics history
        self._add_to_history("neutralization_metrics", metrics)
        
        return metrics
    
    def calculate_pareto_front(self, designs: List[Dict], objectives: List[str]) -> List[int]:
        """
        Calculate the Pareto front for multi-objective optimization.
        
        Args:
            designs: List of design dictionaries
            objectives: List of objective keys to optimize (assumed to be maximize)
            
        Returns:
            List of indices in the Pareto front
        """
        if not designs or not objectives:
            return []
            
        # Extract objective values
        values = np.zeros((len(designs), len(objectives)))
        
        # Handle missing values and objective direction
        for i, design in enumerate(designs):
            for j, obj in enumerate(objectives):
                # Determine if this is a minimize or maximize objective
                minimize = False
                if obj.startswith("minimize_"):
                    obj = obj[9:]  # Remove the minimize_ prefix
                    minimize = True
                    
                # Extract value with nested dictionary support (obj can be "key1.key2")
                value = design
                for key in obj.split("."):
                    if isinstance(value, dict) and key in value:
                        value = value[key]
                    else:
                        value = 0  # Default value if key not found
                        break
                        
                # Convert to float, invert if minimizing
                try:
                    value = float(value)
                    if minimize:
                        value = -value  # Invert for minimization
                except (ValueError, TypeError):
                    value = 0
                    
                values[i, j] = value
        
        # Calculate Pareto front
        pareto_front = []
        for i, design_values in enumerate(values):
            dominated = False
            for j, other_values in enumerate(values):
                if i != j:
                    # Check if design_values is dominated by other_values
                    if np.all(other_values >= design_values) and np.any(other_values > design_values):
                        dominated = True
                        break
            
            if not dominated:
                pareto_front.append(i)
                
        return pareto_front
    
    def track_experiment(self, experiment_name: str, metrics: Dict[str, Any],
                        metadata: Optional[Dict[str, Any]] = None) -> str:
        """
        Track an experiment with associated metrics and metadata.
        
        Args:
            experiment_name: Name of the experiment
            metrics: Dictionary of metrics
            metadata: Optional dictionary of metadata
            
        Returns:
            Experiment ID
        """
        # Create unique experiment ID
        experiment_id = f"{experiment_name}_{self.current_run_id}_{int(time.time())}"  
        
        # Create experiment record
        experiment = {
            "id": experiment_id,
            "name": experiment_name,
            "timestamp": datetime.now().isoformat(),
            "metrics": metrics,
            "metadata": metadata or {}
        }
        
        # Save experiment to disk
        self._save_experiment(experiment)
        
        logger.info(f"Tracked experiment {experiment_id} with {len(metrics)} metrics")
        return experiment_id
    
    def compare_experiments(self, experiment_ids: List[str], 
                         metric_keys: Optional[List[str]] = None) -> Dict[str, Dict[str, List[float]]]:
        """
        Compare metrics across multiple experiments.
        
        Args:
            experiment_ids: List of experiment IDs to compare
            metric_keys: Optional list of specific metric keys to compare
            
        Returns:
            Dictionary mapping metric keys to values across experiments
        """
        comparisons = defaultdict(dict)
        
        # Load each experiment
        for exp_id in experiment_ids:
            experiment = self._load_experiment(exp_id)
            if experiment and "metrics" in experiment:
                for key, value in experiment["metrics"].items():
                    if metric_keys is None or key in metric_keys:
                        if key not in comparisons:
                            comparisons[key] = {}
                        comparisons[key][exp_id] = value
        
        return dict(comparisons)
    
    def visualize_metrics_history(self, metric_type: str, metric_keys: List[str], 
                                output_path: Optional[str] = None) -> str:
        """
        Visualize metrics history over time.
        
        Args:
            metric_type: Type of metrics to visualize (e.g., "binding_metrics")
            metric_keys: List of specific metric keys to visualize
            output_path: Optional output file path
            
        Returns:
            Path to saved visualization or empty string if unsuccessful
        """
        if metric_type not in self.metrics_history or not self.metrics_history[metric_type]:
            logger.warning(f"No history available for metric type: {metric_type}")
            return ""
            
        # Create figure
        plt.figure(figsize=(12, 8))
        
        # Extract timestamps and metrics
        history = self.metrics_history[metric_type]
        timestamps = [entry["timestamp"] for entry in history]
        x_values = range(len(timestamps))  # Use indices for x-axis
        
        # Plot each metric
        for key in metric_keys:
            values = [entry["metrics"].get(key, np.nan) for entry in history]
            plt.plot(x_values, values, marker='o', label=key)
            
        # Set labels and title
        plt.xlabel("Run Index")
        plt.ylabel("Metric Value")
        plt.title(f"{metric_type.replace('_', ' ').title()} History")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Add timestamp ticks (but not too many)
        if len(timestamps) > 1:
            tick_indices = [0, len(timestamps)-1]  # Always show first and last
            if len(timestamps) > 10:
                # Add some intermediate ticks
                step = len(timestamps) // 5
                tick_indices.extend(range(step, len(timestamps), step))
                tick_indices = sorted(set(tick_indices))  # Remove duplicates and sort
            
            plt.xticks(tick_indices, [timestamps[i].split('T')[0] for i in tick_indices], rotation=45)
        
        # Save figure if output path provided
        if output_path:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            plt.tight_layout()
            plt.savefig(output_path)
            plt.close()
            return output_path
        else:
            # Create default output path
            default_path = os.path.join(self.metrics_dir, 
                                     f"{metric_type}_history_{self.current_run_id}.png")
            os.makedirs(os.path.dirname(default_path), exist_ok=True)
            plt.tight_layout()
            plt.savefig(default_path)
            plt.close()
            return default_path
    
    def export_metrics_summary(self, output_path: Optional[str] = None) -> str:
        """
        Export a summary of all tracked metrics.
        
        Args:
            output_path: Optional output file path
            
        Returns:
            Path to saved summary file
        """
        # Create summary dictionary
        summary = {
            "run_id": self.current_run_id,
            "timestamp": datetime.now().isoformat(),
            "metrics_history": dict(self.metrics_history)
        }
        
        # Determine output path
        if not output_path:
            output_path = os.path.join(self.metrics_dir, f"metrics_summary_{self.current_run_id}.json")
            
        # Ensure directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Write summary to file
        try:
            with open(output_path, 'w') as f:
                json.dump(summary, f, indent=2)
            logger.info(f"Exported metrics summary to {output_path}")
            return output_path
        except Exception as e:
            logger.error(f"Failed to export metrics summary: {str(e)}")
            return ""
    
    def _add_to_history(self, metric_type: str, metrics: Dict[str, Any]) -> None:
        """
        Add metrics to the history with timestamp.
        
        Args:
            metric_type: Type of metrics
            metrics: Dictionary of metric values
        """
        entry = {
            "timestamp": datetime.now().isoformat(),
            "metrics": metrics
        }
        
        self.metrics_history[metric_type].append(entry)
    
    def _save_experiment(self, experiment: Dict[str, Any]) -> None:
        """
        Save an experiment record to disk.
        
        Args:
            experiment: Dictionary containing experiment data
        """
        experiment_dir = os.path.join(self.metrics_dir, "experiments")
        os.makedirs(experiment_dir, exist_ok=True)
        
        file_path = os.path.join(experiment_dir, f"{experiment['id']}.json")
        
        try:
            with open(file_path, 'w') as f:
                json.dump(experiment, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save experiment: {str(e)}")
    
    def _load_experiment(self, experiment_id: str) -> Optional[Dict[str, Any]]:
        """
        Load an experiment record from disk.
        
        Args:
            experiment_id: ID of the experiment to load
            
        Returns:
            Dictionary containing experiment data, or None if not found
        """
        file_path = os.path.join(self.metrics_dir, "experiments", f"{experiment_id}.json")
        
        if not os.path.exists(file_path):
            logger.warning(f"Experiment file not found: {file_path}")
            return None
            
        try:
            with open(file_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            logger.error(f"Failed to load experiment: {str(e)}")
            return None
    
    def _calculate_hits_at_k(self, predicted: List[float], true: List[float], k: int) -> float:
        """
        Calculate hits@k metric for binding affinity predictions.
        
        Args:
            predicted: Predicted affinity values
            true: True affinity values
            k: Top-k threshold
            
        Returns:
            Hits@k score (0-1)
        """
        if len(predicted) <= k:
            return 1.0  # All predictions are in top-k
            
        # Create index arrays
        true_sorted_idx = np.argsort(true)
        pred_sorted_idx = np.argsort(predicted)
        
        # Get top-k indices
        true_top_k = set(true_sorted_idx[:k])
        pred_top_k = set(pred_sorted_idx[:k])
        
        # Calculate overlap
        overlap = len(true_top_k.intersection(pred_top_k))
        
        return overlap / k
    
    def _max_sequence_similarity(self, sequence: str, reference_sequences: List[str]) -> float:
        """
        Calculate maximum similarity between a sequence and reference sequences.
        
        Args:
            sequence: Query sequence
            reference_sequences: List of reference sequences
            
        Returns:
            Maximum similarity score (0-1)
        """
        max_similarity = 0.0
        
        for ref_seq in reference_sequences:
            # Simple exact matching for demonstration
            # In practice, would use alignment algorithms
            
            # Ensure sequences are same length for comparison
            min_len = min(len(sequence), len(ref_seq))
            seq = sequence[:min_len]
            ref = ref_seq[:min_len]
            
            # Count matching positions
            matches = sum(a == b for a, b in zip(seq, ref))
            similarity = matches / min_len
            
            max_similarity = max(max_similarity, similarity)
        
        return max_similarity


# Standalone utility functions for validation

def calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Calculate Root Mean Square Deviation (RMSD) between two coordinate sets.
    
    Args:
        coords1: First set of 3D coordinates (N x 3)
        coords2: Second set of 3D coordinates (N x 3)
        
    Returns:
        RMSD value
    """
    if coords1.shape != coords2.shape or coords1.shape[1] != 3:
        raise ValueError("Coordinate arrays must have same shape and be Nx3")
        
    # Calculate squared differences
    squared_diff = np.sum((coords1 - coords2) ** 2, axis=1)
    
    # Calculate RMSD
    rmsd = np.sqrt(np.mean(squared_diff))
    
    return float(rmsd)


def calculate_matthew_correlation(y_true: List[int], y_pred: List[int]) -> float:
    """
    Calculate Matthew's Correlation Coefficient for binary classification.
    
    Args:
        y_true: List of true binary labels (0/1)
        y_pred: List of predicted binary labels (0/1)
        
    Returns:
        MCC value (-1 to 1)
    """
    if len(y_true) != len(y_pred) or len(y_true) == 0:
        logger.error("Labels must be non-empty and of same length")
        return 0.0
    
    # Calculate confusion matrix
    true_pos = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 1)
    true_neg = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 0)
    false_pos = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 1)
    false_neg = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 0)
    
    # Calculate MCC
    numerator = (true_pos * true_neg) - (false_pos * false_neg)
    denominator = ((true_pos + false_pos) * 
                  (true_pos + false_neg) * 
                  (true_neg + false_pos) * 
                  (true_neg + false_neg))
    
    # Handle division by zero
    if denominator == 0:
        return 0.0
        
    return numerator / np.sqrt(denominator)


def evaluate_binding_prediction_model(model_outputs: List[Dict], 
                                   ground_truth: List[Dict]) -> Dict[str, float]:
    """
    Evaluate a binding prediction model against ground truth.
    
    Args:
        model_outputs: List of model output dictionaries
        ground_truth: List of ground truth dictionaries
        
    Returns:
        Dictionary of evaluation metrics
    """
    # Extract predictions and true values
    try:
        predictions = [output.get("binding_score", 0) for output in model_outputs]
        true_values = [truth.get("experimental_binding", 0) for truth in ground_truth]
        
        # For binary classification metrics, convert to binary using threshold
        threshold = 0.5  # Binding threshold
        binary_pred = [1 if p >= threshold else 0 for p in predictions]
        binary_true = [1 if t >= threshold else 0 for t in true_values]
        
        # Calculate metrics
        metrics = {
            "pearson_r": float(stats.pearsonr(predictions, true_values)[0]),
            "spearman_rho": float(stats.spearmanr(predictions, true_values)[0]),
            "mcc": calculate_matthew_correlation(binary_true, binary_pred),
            "accuracy": sum(1 for p, t in zip(binary_pred, binary_true) if p == t) / len(binary_pred),
        }
        
        return metrics
        
    except Exception as e:
        logger.error(f"Error evaluating binding prediction model: {str(e)}")
        return {"error": str(e)}


def load_validation_dataset(file_path: str) -> Tuple[List[Dict], List[Dict]]:
    """
    Load validation dataset from file.
    
    Args:
        file_path: Path to validation dataset file
        
    Returns:
        Tuple of (inputs, ground_truth)
    """
    if not os.path.exists(file_path):
        logger.error(f"Validation dataset not found: {file_path}")
        return [], []
    
    try:
        # Determine file type
        if file_path.endswith(".json"):
            with open(file_path, 'r') as f:
                data = json.load(f)
            return data.get("inputs", []), data.get("ground_truth", [])
            
        elif file_path.endswith(".csv"):
            df = pd.read_csv(file_path)
            
            # Extract inputs and ground truth based on column names
            inputs = []
            ground_truth = []
            
            for _, row in df.iterrows():
                input_dict = {}
                truth_dict = {}
                
                for col in df.columns:
                    if col.startswith("input_"):
                        key = col[6:]  # Remove "input_" prefix
                        input_dict[key] = row[col]
                    elif col.startswith("truth_"):
                        key = col[6:]  # Remove "truth_" prefix
                        truth_dict[key] = row[col]
                
                inputs.append(input_dict)
                ground_truth.append(truth_dict)
                
            return inputs, ground_truth
            
        else:
            logger.error(f"Unsupported file format: {file_path}")
            return [], []
            
    except Exception as e:
        logger.error(f"Error loading validation dataset: {str(e)}")
        return [], []
