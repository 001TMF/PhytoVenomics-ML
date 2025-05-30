#!/usr/bin/env python
# antibody_design/affinity_optimizer.py

import logging
import random
import numpy as np
from typing import Dict, List, Tuple, Optional
from .antibody_generator import HumanAntibodyDesigner

logger = logging.getLogger("Phytovenomics.AffinityOptimizer")

class AffinityOptimizer:
    """
    Optimizes the binding affinity of designed antibodies against toxin epitopes.
    """
    def __init__(self, config: Dict):
        """
        Initialize the AffinityOptimizer with configuration.
        
        Args:
            config: Dictionary containing configuration parameters
        """
        self.config = config
        self.antibody_designer = None  # Will be set when needed
        self.optimization_iterations = self.config.get("optimization_iterations", 3)
        self.optimization_population = self.config.get("optimization_population", 5)
        
    def set_antibody_designer(self, designer: HumanAntibodyDesigner):
        """
        Set the antibody designer instance to use for optimizations.
        
        Args:
            designer: HumanAntibodyDesigner instance
        """
        self.antibody_designer = designer
    
    def optimize_antibody(self, antibody: Dict, target: Dict) -> Dict:
        """
        Optimize the binding affinity of an antibody against a target epitope.
        
        Args:
            antibody: Dictionary containing antibody data
            target: Dictionary containing epitope information
            
        Returns:
            Optimized antibody dictionary
        """
        logger.info(f"Starting affinity optimization for antibody {antibody['id']} against {target['id']}")
        
        # Ensure we have an antibody designer
        if not self.antibody_designer:
            raise ValueError("Antibody designer must be set before optimization")
        
        best_antibody = antibody
        best_energy = antibody.get("predicted_binding", {}).get("energy", 0)
        
        # Run multiple rounds of optimization
        for i in range(self.optimization_iterations):
            logger.info(f"Optimization iteration {i+1}/{self.optimization_iterations}")
            
            # Generate optimized variant
            optimized = self.antibody_designer.cdr_optimization(best_antibody, target)
            
            # Compare with current best
            current_energy = optimized.get("predicted_binding", {}).get("energy", 0)
            
            if current_energy < best_energy:  # Lower energy = better binding
                logger.info(f"Found improved variant: {current_energy} vs {best_energy}")
                best_antibody = optimized
                best_energy = current_energy
            else:
                logger.info("No improvement in this iteration")
        
        # Add optimization metadata
        best_antibody["optimization_info"] = {
            "iterations": self.optimization_iterations,
            "initial_energy": antibody.get("predicted_binding", {}).get("energy", 0),
            "final_energy": best_energy,
            "improvement": antibody.get("predicted_binding", {}).get("energy", 0) - best_energy
        }
        
        return best_antibody
    
    def batch_optimize(self, antibodies: List[Dict], targets: Dict[str, Dict]) -> List[Dict]:
        """
        Optimize a batch of antibodies against their respective targets.
        
        Args:
            antibodies: List of antibody dictionaries
            targets: Dictionary mapping target IDs to target dictionaries
            
        Returns:
            List of optimized antibody dictionaries
        """
        logger.info(f"Batch optimizing {len(antibodies)} antibodies")
        optimized_antibodies = []
        
        for antibody in antibodies:
            target_id = antibody.get("target")
            if target_id in targets:
                optimized = self.optimize_antibody(antibody, targets[target_id])
                optimized_antibodies.append(optimized)
            else:
                logger.warning(f"Target {target_id} not found for antibody {antibody.get('id')}")
                optimized_antibodies.append(antibody)  # Return unoptimized
        
        return optimized_antibodies
    
    def optimize_for_cross_reactivity(self, antibody: Dict, primary_target: Dict, related_targets: List[Dict]) -> Dict:
        """
        Optimize an antibody to bind to multiple related targets (for broadly neutralizing antibodies).
        
        Args:
            antibody: Dictionary containing antibody data
            primary_target: Dictionary containing primary epitope information
            related_targets: List of dictionaries containing information about related epitopes
            
        Returns:
            Optimized antibody dictionary with cross-reactivity information
        """
        logger.info(f"Optimizing antibody {antibody['id']} for cross-reactivity across {len(related_targets)+1} targets")
        
        # First, optimize for the primary target
        optimized_antibody = self.optimize_antibody(antibody, primary_target)
        
        # Calculate binding to related targets
        cross_reactivity = []
        
        for target in related_targets:
            # Extract antibody sequences
            heavy_chain = optimized_antibody["heavy_chain"]["sequence"]
            light_chain = optimized_antibody["light_chain"]["sequence"]
            
            # Predict binding to this target
            binding_energy = self.antibody_designer._predict_binding_energy(
                heavy_chain, light_chain, target["sequence"]
            )
            
            # Calculate affinity
            affinity = self.antibody_designer._convert_energy_to_affinity(binding_energy)
            
            cross_reactivity.append({
                "target_id": target["id"],
                "energy": binding_energy,
                "affinity_nm": affinity
            })
        
        # Add cross-reactivity data to the antibody
        optimized_antibody["cross_reactivity"] = cross_reactivity
        
        # Calculate a cross-reactivity score
        energies = [cr["energy"] for cr in cross_reactivity]
        if energies:
            # Lower mean = better binding across targets
            cross_score = -10 * (1.0 / (1.0 + np.std(energies))) * np.mean(energies)
            optimized_antibody["cross_reactivity_score"] = round(cross_score, 2)
        
        return optimized_antibody
    
    def optimize_cocktail_coverage(self, antibodies: List[Dict], targets: Dict[str, Dict]) -> List[Dict]:
        """
        Optimize a set of antibodies as a cocktail to maximize coverage of multiple targets.
        
        Args:
            antibodies: List of antibody dictionaries
            targets: Dictionary mapping target IDs to target dictionaries
            
        Returns:
            List of optimized antibody dictionaries designed for cocktail usage
        """
        logger.info(f"Optimizing cocktail coverage for {len(antibodies)} antibodies against {len(targets)} targets")
        
        # First, batch optimize all antibodies individually
        optimized_antibodies = self.batch_optimize(antibodies, targets)
        
        # Then analyze for coverage gaps
        coverage = self._analyze_target_coverage(optimized_antibodies, targets)
        
        # For antibodies targeting poorly covered targets, apply more aggressive optimization
        for i, antibody in enumerate(optimized_antibodies):
            target_id = antibody.get("target")
            if target_id in coverage and coverage[target_id] < 0.7:  # Coverage below 70%
                logger.info(f"Applying aggressive optimization for poorly covered target {target_id}")
                
                # Double the optimization iterations for this antibody
                saved_iterations = self.optimization_iterations
                self.optimization_iterations *= 2
                
                # Re-optimize
                optimized_antibodies[i] = self.optimize_antibody(antibody, targets[target_id])
                
                # Restore original iterations
                self.optimization_iterations = saved_iterations
        
        # Add cocktail suitability scores
        for i, antibody in enumerate(optimized_antibodies):
            # Calculate how well this antibody complements the others in the cocktail
            complement_score = self._calculate_complement_score(antibody, optimized_antibodies)
            optimized_antibodies[i]["cocktail_metrics"] = {
                "complement_score": complement_score,
                "coverage_contribution": self._calculate_coverage_contribution(antibody, targets)
            }
        
        return optimized_antibodies
    
    def _analyze_target_coverage(self, antibodies: List[Dict], targets: Dict[str, Dict]) -> Dict[str, float]:
        """Analyze how well the current set of antibodies covers the targets"""
        coverage = {}
        
        for target_id, target in targets.items():
            # Find antibodies targeting this target
            target_antibodies = [ab for ab in antibodies if ab.get("target") == target_id]
            
            if target_antibodies:
                # Calculate average binding energy (lower is better)
                energies = [ab.get("predicted_binding", {}).get("energy", 0) for ab in target_antibodies]
                avg_energy = sum(energies) / len(energies) if energies else 0
                
                # Convert to coverage score (0-1, higher is better)
                coverage[target_id] = max(0, min(1, 0.8 - avg_energy / 10))
            else:
                coverage[target_id] = 0.0
        
        return coverage
    
    def _calculate_complement_score(self, antibody: Dict, all_antibodies: List[Dict]) -> float:
        """Calculate how well this antibody complements others in the cocktail"""
        # This is a simplified implementation
        # In practice, would consider epitope overlap, binding mechanisms, etc.
        
        # Check if this is a broadly neutralizing antibody
        is_bnab = antibody.get("type", "") == "broadly_neutralizing"
        
        # BNAbs generally complement specific antibodies well
        if is_bnab:
            return 0.8 + 0.2 * random.random()  # 0.8-1.0
        
        # For specific antibodies, check how many others target the same epitope
        target_id = antibody.get("target")
        similar_antibodies = [ab for ab in all_antibodies 
                             if ab.get("target") == target_id and ab.get("id") != antibody.get("id")]
        
        # If this is the only antibody targeting this epitope, it's very complementary
        if not similar_antibodies:
            return 0.9 + 0.1 * random.random()  # 0.9-1.0
        
        # Otherwise, it's less complementary (redundant)
        return 0.5 + 0.3 * random.random()  # 0.5-0.8
    
    def _calculate_coverage_contribution(self, antibody: Dict, targets: Dict[str, Dict]) -> float:
        """Calculate how much this antibody contributes to overall target coverage"""
        # Simplified implementation
        # In practice, would analyze epitope overlap, epitope conservation, etc.
        
        target_id = antibody.get("target")
        if target_id not in targets:
            return 0.0
        
        # Get binding energy and specificity
        binding_energy = antibody.get("predicted_binding", {}).get("energy", 0)
        specificity = antibody.get("specificity", 0.5)
        
        # For BNAbs, check cross-reactivity
        if antibody.get("type", "") == "broadly_neutralizing" and "cross_reactivity" in antibody:
            cross_reactive_targets = len(antibody["cross_reactivity"])
            cross_score = antibody.get("cross_reactivity_score", 0) 
            
            # BNAbs with good cross-reactivity contribute more to coverage
            return 0.1 * cross_reactive_targets + 0.5 * (1 - binding_energy / -10) + 0.4 * cross_score
        
        # For specific antibodies, contribution depends on binding strength and specificity
        return 0.6 * (1 - binding_energy / -10) + 0.4 * specificity