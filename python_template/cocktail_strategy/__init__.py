#!/usr/bin/env python
# cocktail_strategy/__init__.py

"""
Package for antibody cocktail optimization algorithms.
This package provides tools for designing optimal combinations of antibodies
to neutralize a set of target toxins.
"""

import os
import logging
import random
import numpy as np
from typing import Dict, List, Any, Union, Optional
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

class CocktailOptimizer:
    """
    Class for optimizing antibody cocktails to neutralize a set of toxins.
    Implements various strategies for selecting antibodies that provide
    maximum coverage of target toxins.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the cocktail optimizer.
        
        Args:
            config: Configuration dictionary with optimization parameters
        """
        self.config = config or {}
        
        # Extract configuration parameters
        self.min_antibodies = self.config.get("min_antibodies", 2)
        self.max_antibodies = self.config.get("max_antibodies", 5)
        self.min_coverage_threshold = self.config.get("min_coverage_threshold", 0.7)
        self.affinity_threshold_nm = self.config.get("affinity_threshold_nm", 50.0)  # nM
        self.prefer_broadly_neutralizing = self.config.get("prefer_broadly_neutralizing", True)
        self.monte_carlo_iterations = self.config.get("monte_carlo_iterations", 1000)
        
        logger.info("Initialized CocktailOptimizer")
    
    def design_optimal_cocktail(self, antibodies: List[Dict], toxins: Dict[str, Dict]) -> Dict[str, Any]:
        """
        Design an optimal antibody cocktail to neutralize a set of toxins.
        
        Args:
            antibodies: List of antibody dictionaries with binding information
            toxins: Dictionary of toxins to target
            
        Returns:
            Dictionary with optimal cocktail information
        """
        logger.info(f"Designing optimal cocktail from {len(antibodies)} antibodies for {len(toxins)} toxins")
        
        # Convert toxins dict to list if needed
        toxin_list = list(toxins.values()) if isinstance(toxins, dict) else toxins
        
        # Analyze the binding profile of each antibody
        binding_matrix = self._create_binding_matrix(antibodies, toxin_list)
        
        # Find optimal cocktail using greedy algorithm
        selected_antibodies, coverage = self._greedy_optimization(antibodies, binding_matrix, toxin_list)
        
        # Calculate average coverage
        coverage_values = list(coverage.values())
        average_coverage = sum(coverage_values) / len(coverage_values) if coverage_values else 0
        
        # Return results
        result = {
            "antibodies": [ab["id"] for ab in selected_antibodies],
            "coverage": coverage,
            "average_coverage": average_coverage,
            "rationale": self._generate_rationale(selected_antibodies, coverage, toxin_list)
        }
        
        logger.info(f"Designed cocktail with {len(selected_antibodies)} antibodies, avg coverage: {average_coverage:.2f}")
        
        return result
    
    def _create_binding_matrix(self, antibodies: List[Dict], toxins: List[Dict]) -> np.ndarray:
        """
        Create a binding matrix between antibodies and toxins.
        
        Args:
            antibodies: List of antibody dictionaries
            toxins: List of toxin dictionaries
            
        Returns:
            2D numpy array of binding strengths (normalized 0-1)
        """
        # Initialize binding matrix
        binding_matrix = np.zeros((len(antibodies), len(toxins)))
        
        # Map toxin IDs to indices
        toxin_ids = [t["id"] if isinstance(t, dict) else t for t in toxins]
        toxin_id_to_index = {tid: i for i, tid in enumerate(toxin_ids)}
        
        # Fill binding matrix
        for ab_idx, antibody in enumerate(antibodies):
            if "binds_to" not in antibody:
                continue
                
            for binding in antibody["binds_to"]:
                toxin_id = binding.get("toxin_id")
                if toxin_id in toxin_id_to_index:
                    # Convert affinity to binding strength (0-1)
                    # Lower affinity (in nM) is better
                    affinity_nm = binding.get("affinity_nm", 1000)
                    
                    # Calculate binding strength using a Hill function-like curve
                    # This gives a value close to 1 for low affinities (strong binding)
                    # and close to 0 for high affinities (weak binding)
                    binding_strength = 1.0 / (1.0 + (affinity_nm / self.affinity_threshold_nm))
                    
                    binding_matrix[ab_idx, toxin_id_to_index[toxin_id]] = binding_strength
        
        return binding_matrix
    
    def _greedy_optimization(self, antibodies: List[Dict], binding_matrix: np.ndarray, 
                            toxins: List[Dict]) -> tuple:
        """
        Use a greedy algorithm to select antibodies for the cocktail.
        
        Args:
            antibodies: List of antibody dictionaries
            binding_matrix: Matrix of binding strengths between antibodies and toxins
            toxins: List of toxin dictionaries
            
        Returns:
            Tuple of (selected antibodies, coverage dictionary)
        """
        # Initialize coverage and selected antibodies
        toxin_ids = [t["id"] if isinstance(t, dict) else t for t in toxins]
        current_coverage = np.zeros(len(toxins))
        selected_antibodies = []
        
        # Adjust antibody weights based on developability and stability
        antibody_weights = np.ones(len(antibodies))
        for i, ab in enumerate(antibodies):
            dev_score = ab.get("developability_score", 0.7)
            stability = ab.get("stability", 0.7)
            broadly_neutralizing = ab.get("is_broadly_neutralizing", False)
            
            # Calculate weight (higher is better)
            weight = 0.4 * dev_score + 0.4 * stability
            
            # Bonus for broadly neutralizing antibodies
            if broadly_neutralizing and self.prefer_broadly_neutralizing:
                weight *= 1.2
                
            antibody_weights[i] = weight
        
        # Iteratively select antibodies greedily
        while len(selected_antibodies) < self.max_antibodies:
            # Calculate the weighted gain for each candidate antibody
            best_gain = 0
            best_antibody_idx = -1
            
            for i in range(len(antibodies)):
                # Skip if already selected
                if antibodies[i] in selected_antibodies:
                    continue
                    
                # Calculate coverage gain with this antibody
                potential_coverage = np.maximum(current_coverage, binding_matrix[i])
                coverage_gain = np.sum(potential_coverage - current_coverage)
                
                # Weight by antibody quality
                weighted_gain = coverage_gain * antibody_weights[i]
                
                if weighted_gain > best_gain:
                    best_gain = weighted_gain
                    best_antibody_idx = i
            
            # Stop if no improvement
            if best_antibody_idx == -1 or best_gain <= 0.01:
                break
                
            # Add the best antibody to the cocktail
            selected_antibodies.append(antibodies[best_antibody_idx])
            current_coverage = np.maximum(current_coverage, binding_matrix[best_antibody_idx])
            
            # Check if we have enough coverage
            if np.min(current_coverage) >= self.min_coverage_threshold and len(selected_antibodies) >= self.min_antibodies:
                break
        
        # Convert coverage to dictionary
        coverage_dict = {toxin_ids[i]: float(current_coverage[i]) for i in range(len(toxins))}
        
        return selected_antibodies, coverage_dict
    
    def _generate_rationale(self, selected_antibodies: List[Dict], 
                           coverage: Dict[str, float], 
                           toxins: List[Dict]) -> str:
        """
        Generate a human-readable rationale for the selected cocktail.
        
        Args:
            selected_antibodies: List of selected antibody dictionaries
            coverage: Dictionary mapping toxin IDs to coverage values
            toxins: List of toxin dictionaries
            
        Returns:
            String with rationale explaining the selection
        """
        # Count broadly neutralizing antibodies
        broadly_neutralizing_count = sum(
            1 for ab in selected_antibodies if ab.get("is_broadly_neutralizing", False)
        )
        
        # Create toxin name mapping
        toxin_names = {}
        for t in toxins:
            if isinstance(t, dict):
                toxin_names[t["id"]] = t.get("name", t["id"])
            else:
                toxin_names[t] = t
        
        # Generate rationale
        lines = [
            f"Selected {len(selected_antibodies)} antibodies including {broadly_neutralizing_count} broadly neutralizing ones.",
            f"Average toxin coverage: {sum(coverage.values()) / len(coverage):.2f}",
            "Antibody selection rationale:"
        ]
        
        for i, ab in enumerate(selected_antibodies):
            ab_id = ab["id"]
            targeting = ab.get("targeting_strategy", "unknown")
            
            # Find toxins this antibody binds to
            bound_toxins = []
            for binding in ab.get("binds_to", []):
                toxin_id = binding.get("toxin_id")
                if toxin_id in coverage:
                    bound_toxins.append(toxin_id)
            
            toxin_str = ", ".join([toxin_names.get(t, t) for t in bound_toxins[:3]])
            if len(bound_toxins) > 3:
                toxin_str += f" and {len(bound_toxins) - 3} more"
                
            lines.append(f"- {ab_id}: {targeting} antibody binding to {toxin_str}")
        
        return "\n".join(lines)
    
    def visualize_cocktail(self, antibodies: List[Dict], 
                         selected_antibodies: List[Dict],
                         toxins: Dict[str, Dict],
                         output_path: Optional[str] = None) -> str:
        """
        Create a visualization of the cocktail performance.
        
        Args:
            antibodies: List of all antibody dictionaries
            selected_antibodies: List of selected antibody dictionaries
            toxins: Dictionary of toxins to target
            output_path: Path to save visualization
            
        Returns:
            Path to saved visualization or empty string if not saved
        """
        # Create binding matrix
        toxin_list = list(toxins.values()) if isinstance(toxins, dict) else toxins
        binding_matrix = self._create_binding_matrix(antibodies, toxin_list)
        
        # Get coverage for selected antibodies
        selected_indices = [antibodies.index(ab) for ab in selected_antibodies]
        coverage = np.max(binding_matrix[selected_indices, :], axis=0)
        
        # Create visualization
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot coverage
        toxin_ids = [t["id"] if isinstance(t, dict) else t for t in toxin_list]
        bars = ax.bar(toxin_ids, coverage, color='skyblue')
        
        # Add threshold line
        ax.axhline(y=self.min_coverage_threshold, linestyle='--', color='red', 
                  label=f'Min threshold ({self.min_coverage_threshold})')
        
        # Add labels and title
        ax.set_xlabel('Toxins')
        ax.set_ylabel('Coverage')
        ax.set_title('Toxin Coverage by Antibody Cocktail')
        ax.set_ylim(0, 1.1)
        
        # Add value labels on top of bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                   f'{height:.2f}', ha='center', va='bottom')
        
        ax.legend()
        plt.tight_layout()
        
        # Save if path provided
        if output_path:
            plt.savefig(output_path)
            plt.close()
            return output_path
        
        return ""

# For backwards compatibility
__version__ = "0.1.0"