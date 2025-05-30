#!/usr/bin/env python
# cocktail_strategy/cocktail_optimizer.py

import logging
import random
import numpy as np
from typing import Dict, List, Optional, Set, Tuple, Any
from collections import defaultdict

logger = logging.getLogger("Phytovenomics.CocktailOptimizer")

class CocktailOptimizer:
    """Optimizer for antibody cocktail composition and strategy."""
    
    def __init__(self, config: Dict):
        """
        Initialize the CocktailOptimizer with configuration.
        
        Args:
            config: Dictionary containing configuration parameters
        """
        self.config = config
        self.max_cocktail_size = config.get("max_cocktail_size", 10)
        self.min_coverage_threshold = config.get("min_coverage_threshold", 0.3)  # Lowered threshold to find optimal cocktails
        self.min_cocktail_size = config.get("min_cocktail_size", 2)  # Ensure at least 2 antibodies
        self.manufacturability_threshold = config.get("manufacturability_threshold", 0.65)
    
    def design_optimal_cocktail(self, antibodies: List[Dict], toxins: Dict[str, Dict]) -> Dict:
        """
        Design an optimal antibody cocktail for neutralizing a set of toxins.
        
        Args:
            antibodies: List of designed antibodies
            toxins: Dictionary mapping toxin IDs to toxin data
            
        Returns:
            Dictionary containing the optimized cocktail design
        """
        logger.info(f"Designing optimal cocktail from {len(antibodies)} antibodies for {len(toxins)} toxins")
        
        # Store reference to all antibodies for later use
        self._all_antibodies = antibodies
        
        # Create binding matrix (antibody x toxin)
        binding_matrix, ab_idx_map, toxin_idx_map = self._create_binding_matrix(antibodies, toxins)
        
        # Group toxins by family
        toxin_families = self._group_toxins_by_family(toxins)
        
        # Select antibody subset using greedy coverage algorithm
        selected_antibodies, coverage_metrics = self._select_antibody_subset(
            binding_matrix, ab_idx_map, toxin_idx_map, toxin_families
        )
        
        # Create cocktail object
        cocktail = {
            "antibodies": selected_antibodies,
            "antibody_count": len(selected_antibodies),
            "coverage_metrics": coverage_metrics,
            "average_coverage": coverage_metrics.get("average_family_coverage", 0),
            "manufacturability_score": self._calculate_manufacturability_score(selected_antibodies),
            "all_antibodies": self._all_antibodies  # Add reference to all antibodies
        }
        
        logger.info(f"Designed cocktail with {len(selected_antibodies)} antibodies, "
                  f"{cocktail['average_coverage']:.2%} average family coverage")
        
        return cocktail
    
    def generate_cocktail_variants(self, base_cocktail: Dict, antibodies: List[Dict], 
                                  toxins: Dict[str, Dict], variant_count: int) -> List[Dict]:
        """
        Generate variants of a cocktail for comparison.
        
        Args:
            base_cocktail: Base cocktail design
            antibodies: List of all available antibodies
            toxins: Dictionary mapping toxin IDs to toxin data
            variant_count: Number of variants to generate
            
        Returns:
            List of cocktail variant dictionaries
        """
        logger.info(f"Generating {variant_count} cocktail variants")
        
        variants = [base_cocktail]  # Start with the base cocktail
        
        # Categorize antibodies
        bnabs = [ab for ab in antibodies if ab.get("is_broadly_neutralizing", False)]
        specific_abs = [ab for ab in antibodies if not ab.get("is_broadly_neutralizing", False)]
        
        # Create binding matrix for scoring
        binding_matrix, ab_idx_map, toxin_idx_map = self._create_binding_matrix(antibodies, toxins)
        
        # Group toxins by family
        toxin_families = self._group_toxins_by_family(toxins)
        
        # Generate variants with different strategies
        strategies = [
            "bnab_heavy",        # Prefer broadly neutralizing antibodies
            "specific_focused",  # Focus on specific, high-affinity antibodies
            "family_coverage",   # Maximize coverage of toxin families
            "mixed_approach"     # Balance of specific and broad antibodies
        ]
        
        # Generate variants with different strategies
        for strategy in strategies[:variant_count]:
            variant = self._generate_variant_by_strategy(
                strategy, antibodies, bnabs, specific_abs, 
                binding_matrix, ab_idx_map, toxin_idx_map, toxin_families
            )
            variants.append(variant)
        
        return variants[:variant_count + 1]  # Include base cocktail + up to variant_count variants
    
    def optimize_for_manufacturability(self, cocktail: Dict) -> Dict:
        """
        Optimize a cocktail design for manufacturability.
        
        Args:
            cocktail: Original cocktail design
            
        Returns:
            Optimized cocktail dictionary
        """
        logger.info("Optimizing cocktail for manufacturability")
        
        # Create a copy of the cocktail to modify
        optimized_cocktail = cocktail.copy()
        selected_antibodies = cocktail.get("antibodies", [])
        
        # Current manufacturability score
        current_score = self._calculate_manufacturability_score(selected_antibodies)
        
        if current_score >= self.manufacturability_threshold:
            logger.info(f"Cocktail already meets manufacturability threshold ({current_score:.2f})")
            optimized_cocktail["manufacturability_score"] = current_score
            return optimized_cocktail
        
        # Identify problematic antibodies (those with poor developability)
        problematic_abs = []
        for ab in selected_antibodies:
            developability = ab.get("developability_score", 0.5)
            if developability < 0.6:  # Threshold for developability concerns
                problematic_abs.append(ab)
        
        if not problematic_abs:
            logger.info("No problematic antibodies identified")
            optimized_cocktail["manufacturability_score"] = current_score
            return optimized_cocktail
        
        # Sort by developability (ascending) - worst first
        problematic_abs.sort(key=lambda x: x.get("developability_score", 0))
        
        # Attempt to replace problematic antibodies with similar ones that have better developability
        # This is simplified - a real implementation would use more sophisticated replacement strategies
        for prob_ab in problematic_abs[:2]:  # Try to replace up to 2 problematic antibodies
            # Find replacement candidates that target the same epitope/toxin
            target_toxin = prob_ab.get("target_toxin_id")
            target_epitope = prob_ab.get("target_epitope_id")
            
            # Check if we have the full antibody pool
            # If not, we'll just adjust the current antibodies
            if "all_antibodies" in optimized_cocktail:
                all_abs = optimized_cocktail["all_antibodies"]
                
                # Find potential replacements
                replacements = [
                    ab for ab in all_abs 
                    if ab.get("target_toxin_id") == target_toxin and
                    ab.get("developability_score", 0) > prob_ab.get("developability_score", 0) and
                    ab.get("id") != prob_ab.get("id")
                ]
                
                if replacements:
                    # Sort by developability (descending)
                    replacements.sort(key=lambda x: x.get("developability_score", 0), reverse=True)
                    
                    # Replace the problematic antibody
                    selected_antibodies = [ab for ab in selected_antibodies if ab.get("id") != prob_ab.get("id")]
                    selected_antibodies.append(replacements[0])
            
        # If we can't replace, try to optimize the existing antibodies
        # In a real implementation, would make targeted mutations to improve developability
        # while maintaining binding properties
        for i, ab in enumerate(selected_antibodies):
            if ab.get("developability_score", 0) < 0.6:
                # Simulate improving developability through protein engineering
                # In reality, this would be a complex process involving mutation and testing
                improved_ab = ab.copy()
                improved_ab["developability_score"] = min(1.0, ab.get("developability_score", 0.5) + 0.15)
                improved_ab["notes"] = "Engineered for improved developability"
                selected_antibodies[i] = improved_ab
        
        # Recalculate manufacturability score
        new_score = self._calculate_manufacturability_score(selected_antibodies)
        
        optimized_cocktail["antibodies"] = selected_antibodies
        optimized_cocktail["manufacturability_score"] = new_score
        optimized_cocktail["notes"] = "Optimized for manufacturability"
        
        logger.info(f"Improved manufacturability score from {current_score:.2f} to {new_score:.2f}")
        
        return optimized_cocktail
    
    def _create_binding_matrix(self, antibodies: List[Dict], toxins: Dict[str, Dict]) -> Tuple[np.ndarray, Dict, Dict]:
        """Create a binding affinity matrix between antibodies and toxins."""
        # Create maps from IDs to indices
        ab_idx_map = {ab["id"]: i for i, ab in enumerate(antibodies)}
        toxin_idx_map = {toxin_id: i for i, toxin_id in enumerate(toxins.keys())}
        
        # Create matrix (rows=antibodies, columns=toxins)
        binding_matrix = np.zeros((len(antibodies), len(toxins)))
        
        # Fill matrix with binding affinities
        for ab_idx, antibody in enumerate(antibodies):
            # Get list of toxins this antibody binds to
            binds_to = antibody.get("binds_to", [])
            
            for binding in binds_to:
                toxin_id = binding.get("toxin_id")
                if toxin_id in toxin_idx_map:
                    toxin_idx = toxin_idx_map[toxin_id]
                    
                    # Get affinity (lower value = stronger binding)
                    # Convert from nM to a 0-1 score (1 = strongest binding)
                    affinity_nm = binding.get("affinity_nm", 1000)
                    binding_score = 1.0 / (1.0 + (affinity_nm / 100.0))  # Sigmoid-like normalization
                    
                    # Ensure binding score is in range [0,1]
                    binding_score = min(1.0, max(0.0, binding_score))
                    
                    binding_matrix[ab_idx, toxin_idx] = binding_score
        
        return binding_matrix, ab_idx_map, toxin_idx_map
    
    def _group_toxins_by_family(self, toxins: Dict[str, Dict]) -> Dict[str, List[str]]:
        """Group toxins by their family."""
        families = defaultdict(list)
        
        for toxin_id, toxin_data in toxins.items():
            family = toxin_data.get("family", "unknown")
            families[family].append(toxin_id)
        
        return dict(families)
    
    def _select_antibody_subset(self, binding_matrix: np.ndarray, ab_idx_map: Dict[str, int], 
                               toxin_idx_map: Dict[str, int], 
                               toxin_families: Dict[str, List[str]]) -> Tuple[List[Dict], Dict]:
        """Select a subset of antibodies to create an optimal cocktail."""
        # Create reverse mappings
        idx_to_ab = {idx: ab_id for ab_id, idx in ab_idx_map.items()}
        idx_to_toxin = {idx: toxin_id for toxin_id, idx in toxin_idx_map.items()}
        
        # Map antibody indices and ensure we handle numeric indices correctly
        antibodies_by_idx = {}
        for ab_id, idx in ab_idx_map.items():
            antibodies_by_idx[idx] = ab_id
        
        # Debug the mappings to find issues
        logger.info(f"Antibody ID to index map: {dict(list(ab_idx_map.items())[:5])}...")
        logger.info(f"Index to antibody ID map: {dict(list(idx_to_ab.items())[:5])}...")
        
        # Create family coverage matrix
        family_coverage_matrix = np.zeros((binding_matrix.shape[0], len(toxin_families)))
        
        for fam_idx, (family, toxin_ids) in enumerate(toxin_families.items()):
            for toxin_id in toxin_ids:
                if toxin_id in toxin_idx_map:
                    toxin_idx = toxin_idx_map[toxin_id]
                    family_coverage_matrix[:, fam_idx] += binding_matrix[:, toxin_idx]
            
            # Normalize by number of toxins in family
            family_coverage_matrix[:, fam_idx] /= len(toxin_ids)
        
        # Greedy selection algorithm
        selected_indices = []
        # Only use indices that are actually in our mapping
        remaining_indices = [i for i in range(binding_matrix.shape[0]) if i in idx_to_ab]
        current_coverage = np.zeros(len(toxin_families))
        
        while len(selected_indices) < self.max_cocktail_size:
            # Calculate incremental coverage for each remaining antibody
            incremental_coverage = np.zeros(len(remaining_indices))
            
            for i, ab_idx in enumerate(remaining_indices):
                # Calculate new coverage if this antibody is added
                new_coverage = np.maximum(current_coverage, family_coverage_matrix[ab_idx])
                
                # Calculate incremental benefit
                incremental_coverage[i] = np.sum(new_coverage - current_coverage)
            
            # If no incremental benefit and we've reached minimum size, we can break
            # Otherwise, continue to add for diversity even with minimal benefit
            if np.max(incremental_coverage) <= 0.01 and len(selected_indices) >= self.min_cocktail_size:  
                break
                
            # Select the antibody with the highest incremental coverage
            if len(incremental_coverage) > 0:  # Only if we have antibodies to choose from
                best_idx = remaining_indices[np.argmax(incremental_coverage)]
                selected_indices.append(best_idx)
                remaining_indices.remove(best_idx)
                
                # Update current coverage
                current_coverage = np.maximum(current_coverage, family_coverage_matrix[best_idx])
            else:
                # No more antibodies to add
                break
            
            # Continue adding antibodies until we meet both thresholds:
            # 1. Have at least min_cocktail_size antibodies
            # 2. Have reached min_coverage_threshold
            if len(selected_indices) >= self.min_cocktail_size and np.mean(current_coverage) >= self.min_coverage_threshold:
                logger.info(f"Sufficient coverage ({np.mean(current_coverage):.2%}) reached with {len(selected_indices)} antibodies")
                break
                
            # Keep adding if we haven't reached min_cocktail_size regardless of coverage
            if len(selected_indices) < self.min_cocktail_size and len(remaining_indices) == 0:
                logger.warning(f"Unable to reach minimum cocktail size of {self.min_cocktail_size}. Selected {len(selected_indices)} antibodies.")
                break
        
        # Create list of selected antibodies
        selected_ab_ids = []
        for idx in selected_indices:
            if idx in idx_to_ab:
                selected_ab_ids.append(idx_to_ab[idx])
            else:
                logger.warning(f"Index {idx} not found in idx_to_ab mapping!")
        
        logger.info(f"Selected indices: {selected_indices}")
        logger.info(f"Corresponding antibody IDs: {selected_ab_ids}")
        
        # Calculate coverage metrics
        family_coverage = {}
        for fam_idx, family in enumerate(toxin_families.keys()):
            family_coverage[family] = float(current_coverage[fam_idx])
        
        coverage_metrics = {
            "family_coverage": family_coverage,
            "average_family_coverage": float(np.mean(list(family_coverage.values()))),
            "min_family_coverage": float(np.min(list(family_coverage.values()))),
            "max_family_coverage": float(np.max(list(family_coverage.values())))
        }
        
        # Debug log for selected indices and antibody IDs
        logger.info(f"Selected {len(selected_indices)} antibodies with indices: {selected_indices}")
        logger.info(f"Selected antibody IDs: {selected_ab_ids}")
        
        # Return selected antibody objects
        selected_antibodies = []
        for ab_id in selected_ab_ids:
            for ab in self._all_antibodies:
                if ab.get("id") == ab_id:
                    selected_antibodies.append(ab)
                    break
        
        return selected_antibodies, coverage_metrics
    
    def _calculate_manufacturability_score(self, antibodies: List[Dict]) -> float:
        """Calculate the manufacturability score of a set of antibodies."""
        if not antibodies:
            return 0.0
            
        # Average developability score
        avg_developability = np.mean([ab.get("developability_score", 0.5) for ab in antibodies])
        
        # Sequence diversity (lower is better for manufacturing)
        seq_diversity_penalty = 0.0
        if len(antibodies) > 1:
            # Calculate pairwise sequence similarity
            similarities = []
            for i in range(len(antibodies)):
                for j in range(i+1, len(antibodies)):
                    ab1_seq = antibodies[i].get("sequence", "")
                    ab2_seq = antibodies[j].get("sequence", "")
                    
                    # Skip if sequences not available
                    if not ab1_seq or not ab2_seq:
                        continue
                        
                    # Calculate similarity (simplified)
                    similarity = self._sequence_similarity(ab1_seq, ab2_seq)
                    similarities.append(similarity)
            
            # Average similarity (higher is better for manufacturing)
            if similarities:
                avg_similarity = np.mean(similarities)
                # Convert to penalty (0-0.2)
                seq_diversity_penalty = 0.2 * (1.0 - avg_similarity)
        
        # CDR length penalty (extreme lengths can cause manufacturing issues)
        cdr_length_penalty = 0.0
        for ab in antibodies:
            cdr_h3_length = ab.get("cdr_lengths", {}).get("CDR-H3", 15) if ab.get("cdr_lengths") else 15
            
            # Penalize very short or very long CDR-H3
            if cdr_h3_length < 8 or cdr_h3_length > 22:
                cdr_length_penalty += 0.05
        
        # Normalize cdr length penalty
        cdr_length_penalty = min(0.2, cdr_length_penalty)
        
        # Scale count penalty (more antibodies = more complex manufacturing)
        count_penalty = 0.1 * min(1.0, (len(antibodies) - 3) / 7) if len(antibodies) > 3 else 0.0
        
        # Combine factors
        manufacturability_score = 0.6 * avg_developability - seq_diversity_penalty - cdr_length_penalty - count_penalty
        
        # Ensure score is in range [0,1]
        return min(1.0, max(0.0, manufacturability_score))
    
    def _sequence_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence similarity between two antibody sequences."""
        if not seq1 or not seq2:
            return 0.0
            
        # Simple exact matching for demonstration
        # In practice, would use alignment algorithms
        
        # Ensure sequences are same length for comparison
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]
        
        # Count matching positions
        matches = sum(a == b for a, b in zip(seq1, seq2))
        
        return matches / min_len
    
    def _generate_variant_by_strategy(self, strategy: str, antibodies: List[Dict], 
                                     bnabs: List[Dict], specific_abs: List[Dict],
                                     binding_matrix: np.ndarray, ab_idx_map: Dict[str, int],
                                     toxin_idx_map: Dict[str, int], 
                                     toxin_families: Dict[str, List[str]]) -> Dict:
        """Generate a cocktail variant using a specific strategy."""
        # Start with a subset based on strategy
        initial_subset = []
        
        if strategy == "bnab_heavy":
            # Start with top broadly neutralizing antibodies
            sorted_bnabs = sorted(bnabs, 
                                  key=lambda x: x.get("breadth_score", 0) * x.get("developability_score", 0.5),
                                  reverse=True)
            initial_subset = sorted_bnabs[:min(3, len(sorted_bnabs))]
            
        elif strategy == "specific_focused":
            # Start with highest affinity specific antibodies
            sorted_specific = sorted(specific_abs,
                                    key=lambda x: x.get("avg_affinity", 1000),
                                    reverse=False)  # Lower affinity value = better
            
            # Select one antibody per toxin family (if available)
            family_covered = set()
            for ab in sorted_specific:
                family = ab.get("target_toxin_family", "")
                if family and family not in family_covered:
                    initial_subset.append(ab)
                    family_covered.add(family)
                    
                if len(initial_subset) >= 4:  # Start with up to 4 specific antibodies
                    break
                    
        elif strategy == "family_coverage":
            # Select one antibody for each toxin family, prioritizing broadly neutralizing
            for family in toxin_families:
                # Find best antibody for this family
                best_ab = None
                best_score = 0
                
                for ab in antibodies:
                    binds_to_family = False
                    ab_family_coverage = 0
                    
                    for binding in ab.get("binds_to", []):
                        toxin_id = binding.get("toxin_id")
                        if toxin_id in toxins and toxins[toxin_id].get("family") == family:
                            binds_to_family = True
                            ab_family_coverage += 1
                    
                    if binds_to_family:
                        # Score based on coverage and developability
                        score = ab_family_coverage * ab.get("developability_score", 0.5)
                        if score > best_score:
                            best_ab = ab
                            best_score = score
                
                if best_ab:
                    initial_subset.append(best_ab)
            
        elif strategy == "mixed_approach":
            # Mix of broadly neutralizing and specific antibodies
            # Add 1-2 top broadly neutralizing antibodies
            if bnabs:
                sorted_bnabs = sorted(bnabs, key=lambda x: x.get("breadth_score", 0), reverse=True)
                initial_subset.extend(sorted_bnabs[:min(2, len(sorted_bnabs))])
            
            # Add 2-3 specific antibodies for uncovered families
            covered_families = set()
            for ab in initial_subset:
                for binding in ab.get("binds_to", []):
                    toxin_id = binding.get("toxin_id")
                    if toxin_id in toxins:
                        covered_families.add(toxins[toxin_id].get("family", ""))
            
            # Add specific antibodies for uncovered families
            for ab in specific_abs:
                if len(initial_subset) >= 5:  # Limit initial subset size
                    break
                    
                ab_family = ""
                for binding in ab.get("binds_to", []):
                    toxin_id = binding.get("toxin_id")
                    if toxin_id in toxins:
                        ab_family = toxins[toxin_id].get("family", "")
                        break
                
                if ab_family and ab_family not in covered_families:
                    initial_subset.append(ab)
                    covered_families.add(ab_family)
        
        # Initialize with the strategy-specific subset
        selected_ab_indices = [ab_idx_map.get(ab.get("id")) for ab in initial_subset if ab.get("id") in ab_idx_map]
        selected_ab_indices = [idx for idx in selected_ab_indices if idx is not None]  # Filter out None values
        
        # Calculate current coverage
        current_coverage = np.zeros(len(toxin_families))
        family_indices = {family: i for i, family in enumerate(toxin_families.keys())}
        
        # Update current coverage with initial subset
        for ab_idx in selected_ab_indices:
            for toxin_idx in range(binding_matrix.shape[1]):
                if binding_matrix[ab_idx, toxin_idx] > 0:
                    toxin_id = list(toxin_idx_map.keys())[list(toxin_idx_map.values()).index(toxin_idx)]
                    if toxin_id in toxins:
                        family = toxins[toxin_id].get("family", "")
                        if family in family_indices:
                            family_idx = family_indices[family]
                            current_coverage[family_idx] = max(current_coverage[family_idx], binding_matrix[ab_idx, toxin_idx])
        
        # Complete the cocktail using greedy algorithm
        remaining_ab_indices = [i for i in range(binding_matrix.shape[0]) if i not in selected_ab_indices]
        
        while len(selected_ab_indices) < self.max_cocktail_size:
            # Calculate incremental coverage for each remaining antibody
            incremental_coverage = np.zeros(len(remaining_ab_indices))
            
            for i, ab_idx in enumerate(remaining_ab_indices):
                # Calculate new coverage if this antibody is added
                new_coverage = current_coverage.copy()
                
                for toxin_idx in range(binding_matrix.shape[1]):
                    if binding_matrix[ab_idx, toxin_idx] > 0:
                        toxin_id = list(toxin_idx_map.keys())[list(toxin_idx_map.values()).index(toxin_idx)]
                        if toxin_id in toxins:
                            family = toxins[toxin_id].get("family", "")
                            if family in family_indices:
                                family_idx = family_indices[family]
                                new_coverage[family_idx] = max(new_coverage[family_idx], binding_matrix[ab_idx, toxin_idx])
                
                # Calculate incremental benefit
                incremental_coverage[i] = np.sum(new_coverage - current_coverage)
            
            # If no incremental benefit, break
            if len(incremental_coverage) == 0 or np.max(incremental_coverage) <= 0.01:
                break
                
            # Select the antibody with the highest incremental coverage
            best_idx_pos = np.argmax(incremental_coverage)
            best_idx = remaining_ab_indices[best_idx_pos]
            selected_ab_indices.append(best_idx)
            remaining_ab_indices.remove(best_idx)
            
            # Update current coverage
            for toxin_idx in range(binding_matrix.shape[1]):
                if binding_matrix[best_idx, toxin_idx] > 0:
                    toxin_id = list(toxin_idx_map.keys())[list(toxin_idx_map.values()).index(toxin_idx)]
                    if toxin_id in toxins:
                        family = toxins[toxin_id].get("family", "")
                        if family in family_indices:
                            family_idx = family_indices[family]
                            current_coverage[family_idx] = max(current_coverage[family_idx], binding_matrix[best_idx, toxin_idx])
            
            # If we've reached sufficient coverage, we can stop
            if np.mean(current_coverage) >= self.min_coverage_threshold and len(selected_ab_indices) >= 3:
                break
        
        # Create list of selected antibodies
        idx_to_ab = {idx: ab_id for ab_id, idx in ab_idx_map.items()}
        selected_ab_ids = [idx_to_ab[idx] for idx in selected_ab_indices]
        
        # Calculate coverage metrics
        family_coverage = {}
        for family, fam_idx in family_indices.items():
            family_coverage[family] = float(current_coverage[fam_idx])
        
        coverage_metrics = {
            "family_coverage": family_coverage,
            "average_family_coverage": float(np.mean(list(family_coverage.values()))),
            "min_family_coverage": float(np.min(list(family_coverage.values()))),
            "max_family_coverage": float(np.max(list(family_coverage.values())))
        }
        
        # Create selected antibody list
        selected_antibodies = []
        for ab_id in selected_ab_ids:
            for ab in antibodies:
                if ab.get("id") == ab_id:
                    selected_antibodies.append(ab)
                    break
        
        # Calculate manufacturability score
        manufacturability_score = self._calculate_manufacturability_score(selected_antibodies)
        
        # Create the cocktail object
        cocktail = {
            "antibodies": selected_antibodies,
            "antibody_count": len(selected_antibodies),
            "strategy": strategy,
            "coverage_metrics": coverage_metrics,
            "average_coverage": coverage_metrics.get("average_family_coverage", 0),
            "manufacturability_score": manufacturability_score,
            "all_antibodies": antibodies  # Keep reference to all antibodies for optimization
        }
        
        return cocktail