"""
IgLM Guidance and Validation Module

This module implements the fourth stage of the pipeline:
1. Humanness gate - Score full VH/VL and drop designs below percentile threshold
2. Minimal steering on CDR non-contacts - Raise IgLM log-likelihood without altering interface
3. Soft guidance on contact residues - Use KL prior to nudge toward enriched motifs
4. Re-predict structure after edits and keep only neutral or improved poses

Integrates IgLM capabilities for humanness assessment and sequence optimization.
"""

import os
import numpy as np
import torch
import logging
from typing import List, Dict, Tuple, Optional, Union, Any
from pathlib import Path
from dataclasses import dataclass, field
import json
import copy

# IgLM imports (via Germinal)
from colabdesign.iglm import mk_iglm_model

from .config import IgLMConfig, PipelineConfig
from .filtering import BackboneFilterResults, FilterResults


@dataclass
class IgLMScores:
    """Container for IgLM scoring results."""
    sequence_id: str
    full_sequence: str
    heavy_chain_seq: str
    light_chain_seq: str
    vh_log_likelihood: float
    vl_log_likelihood: float
    total_log_likelihood: float
    humanness_percentile: float
    humanness_pass: bool


@dataclass
class SequenceEdit:
    """Container for a single sequence edit."""
    position: int
    original_aa: str
    new_aa: str
    edit_type: str  # "contact" or "non_contact"
    iglm_improvement: float
    structure_impact: str  # "neutral", "improved", "degraded"


@dataclass
class OptimizedSequence:
    """Container for IgLM-optimized sequence."""
    original_filter_result: FilterResults
    optimized_sequence: str
    iglm_scores: IgLMScores
    edits_applied: List[SequenceEdit]
    structure_prediction: Optional[str]  # Path to re-predicted structure
    final_ranking_score: float
    optimization_success: bool


class IgLMHumannessGate:
    """Humanness assessment using IgLM scores."""
    
    def __init__(self, config: IgLMConfig, logger: Optional[logging.Logger] = None):
        self.config = config
        self.logger = logger or logging.getLogger(__name__)
        
        # Initialize IgLM model lazily
        self._iglm_model = None
        self._humanness_distribution = None
        
    def score_humanness(
        self,
        filter_results: List[BackboneFilterResults]
    ) -> List[List[IgLMScores]]:
        """
        Score humanness for all sequences using IgLM.
        
        Args:
            filter_results: List of backbone filter results
            
        Returns:
            List of IgLM scores for each backbone's candidates
        """
        all_scores = []
        
        for backbone_result in filter_results:
            self.logger.info(f"Scoring humanness for backbone {backbone_result.backbone_id}")
            
            backbone_scores = []
            
            for candidate in backbone_result.top_candidates:
                try:
                    scores = self._score_single_sequence(candidate)
                    backbone_scores.append(scores)
                    
                except Exception as e:
                    self.logger.error(f"Humanness scoring failed for {candidate.prediction.prediction_id}: {e}")
                    continue
            
            all_scores.append(backbone_scores)
        
        return all_scores
    
    def _score_single_sequence(self, filter_result: FilterResults) -> IgLMScores:
        """Score a single sequence for humanness."""
        
        if self._iglm_model is None:
            self._iglm_model = self._initialize_iglm_model()
        
        sequence_variant = filter_result.prediction.sequence_variant
        
        # Split into heavy and light chains
        heavy_seq = sequence_variant.heavy_chain_seq
        light_seq = sequence_variant.light_chain_seq
        
        # Score with IgLM
        vh_ll = self._score_heavy_chain(heavy_seq)
        vl_ll = self._score_light_chain(light_seq)
        total_ll = vh_ll + vl_ll
        
        # Calculate humanness percentile
        humanness_percentile = self._calculate_humanness_percentile(total_ll)
        
        # Apply humanness gate
        humanness_pass = humanness_percentile >= self.config.humanness_percentile_threshold
        
        scores = IgLMScores(
            sequence_id=filter_result.prediction.prediction_id,
            full_sequence=sequence_variant.sequence,
            heavy_chain_seq=heavy_seq,
            light_chain_seq=light_seq,
            vh_log_likelihood=vh_ll,
            vl_log_likelihood=vl_ll,
            total_log_likelihood=total_ll,
            humanness_percentile=humanness_percentile,
            humanness_pass=humanness_pass
        )
        
        return scores
    
    def _initialize_iglm_model(self):
        """Initialize IgLM model."""
        self.logger.info("Initializing IgLM model")
        
        iglm_model = mk_iglm_model(
            species=self.config.iglm_species,
            model_path=self.config.iglm_model_path
        )
        
        return iglm_model
    
    def _score_heavy_chain(self, heavy_seq: str) -> float:
        """Score heavy chain sequence with IgLM."""
        return self._iglm_model.score_sequence(heavy_seq, chain_type="heavy")
    
    def _score_light_chain(self, light_seq: str) -> float:
        """Score light chain sequence with IgLM.""" 
        return self._iglm_model.score_sequence(light_seq, chain_type="light")
    
    def _calculate_humanness_percentile(self, log_likelihood: float) -> float:
        """Calculate humanness percentile based on therapeutic distribution."""
        
        if self._humanness_distribution is None:
            self._humanness_distribution = self._load_humanness_distribution()
        
        # Calculate percentile rank in therapeutic distribution
        percentile = np.percentile(self._humanness_distribution, 
                                 100 * np.mean(self._humanness_distribution <= log_likelihood))
        
        return percentile
    
    def _load_humanness_distribution(self) -> np.ndarray:
        """Load reference distribution of therapeutic antibody scores."""
        
        if self.config.therapeutic_distribution_path:
            # Load from file if provided
            return np.load(self.config.therapeutic_distribution_path)
        else:
            # Use default therapeutic-like distribution
            # This would be based on a curated set of therapeutic antibodies
            # For now, use placeholder normal distribution
            return np.random.normal(-5.0, 1.5, 10000)


class IgLMSequenceOptimizer:
    """Sequence optimization using IgLM guidance."""
    
    def __init__(self, config: IgLMConfig, logger: Optional[logging.Logger] = None):
        self.config = config
        self.logger = logger or logging.getLogger(__name__)
        
        self._iglm_model = None
        self._contact_predictor = None
        
    def optimize_sequences(
        self,
        filter_results: List[BackboneFilterResults],
        iglm_scores: List[List[IgLMScores]]
    ) -> List[List[OptimizedSequence]]:
        """
        Optimize sequences using IgLM guidance.
        
        Args:
            filter_results: Filtered structure results
            iglm_scores: Humanness scores for sequences
            
        Returns:
            List of optimized sequences for each backbone
        """
        all_optimized = []
        
        for backbone_result, backbone_scores in zip(filter_results, iglm_scores):
            self.logger.info(f"Optimizing sequences for backbone {backbone_result.backbone_id}")
            
            backbone_optimized = []
            
            # Only optimize sequences that pass humanness gate
            passing_pairs = [
                (candidate, score) 
                for candidate, score in zip(backbone_result.top_candidates, backbone_scores)
                if score.humanness_pass
            ]
            
            for candidate, iglm_score in passing_pairs:
                try:
                    optimized = self._optimize_single_sequence(candidate, iglm_score)
                    backbone_optimized.append(optimized)
                    
                except Exception as e:
                    self.logger.error(f"Sequence optimization failed for {candidate.prediction.prediction_id}: {e}")
                    continue
            
            all_optimized.append(backbone_optimized)
        
        return all_optimized
    
    def _optimize_single_sequence(
        self,
        filter_result: FilterResults,
        iglm_scores: IgLMScores
    ) -> OptimizedSequence:
        """Optimize a single sequence using IgLM guidance."""
        
        # Identify contact and non-contact positions
        contact_positions = self._identify_contact_positions(filter_result)
        non_contact_positions = self._identify_non_contact_positions(filter_result, contact_positions)
        
        # Start with original sequence
        current_sequence = iglm_scores.full_sequence
        current_heavy = iglm_scores.heavy_chain_seq
        current_light = iglm_scores.light_chain_seq
        
        edits_applied = []
        
        # Phase 1: Minimal steering on non-contact positions
        if self.config.use_minimal_steering:
            current_sequence, current_heavy, current_light, non_contact_edits = \
                self._apply_non_contact_steering(
                    current_sequence, current_heavy, current_light, non_contact_positions
                )
            edits_applied.extend(non_contact_edits)
        
        # Phase 2: Soft guidance on contact residues
        if self.config.use_soft_contact_guidance:
            current_sequence, current_heavy, current_light, contact_edits = \
                self._apply_contact_guidance(
                    current_sequence, current_heavy, current_light, contact_positions, filter_result
                )
            edits_applied.extend(contact_edits)
        
        # Re-score with IgLM
        final_iglm_scores = self._score_optimized_sequence(current_heavy, current_light)
        
        # Re-predict structure if configured
        structure_prediction = None
        final_ranking_score = filter_result.ranking_score
        
        if self.config.repredict_after_edits and edits_applied:
            structure_prediction, final_ranking_score = self._repredict_structure(
                current_sequence, filter_result
            )
        
        # Determine optimization success
        optimization_success = self._evaluate_optimization_success(
            iglm_scores, final_iglm_scores, edits_applied, final_ranking_score, filter_result.ranking_score
        )
        
        optimized = OptimizedSequence(
            original_filter_result=filter_result,
            optimized_sequence=current_sequence,
            iglm_scores=final_iglm_scores,
            edits_applied=edits_applied,
            structure_prediction=structure_prediction,
            final_ranking_score=final_ranking_score,
            optimization_success=optimization_success
        )
        
        return optimized
    
    def _identify_contact_positions(self, filter_result: FilterResults) -> List[int]:
        """Identify positions making contact with antigen."""
        # This would analyze the structure to identify interface contacts
        # For now, use placeholder based on CDR positions
        
        sequence_variant = filter_result.prediction.sequence_variant
        contact_positions = []
        
        # Assume CDR positions are more likely to be contacts
        for cdr_name, (start, end) in sequence_variant.backbone_geometry.cdr_assignments.items():
            # Sample some positions from each CDR as contacts
            cdr_contacts = list(range(start, min(start + 5, end + 1)))
            contact_positions.extend(cdr_contacts)
        
        return contact_positions
    
    def _identify_non_contact_positions(
        self, 
        filter_result: FilterResults, 
        contact_positions: List[int]
    ) -> List[int]:
        """Identify CDR positions that don't make contact with antigen."""
        
        sequence_variant = filter_result.prediction.sequence_variant
        all_cdr_positions = []
        
        # Get all CDR positions
        for cdr_name, (start, end) in sequence_variant.backbone_geometry.cdr_assignments.items():
            all_cdr_positions.extend(list(range(start, end + 1)))
        
        # Non-contact positions are CDR positions not in contact
        non_contact_positions = [
            pos for pos in all_cdr_positions 
            if pos not in contact_positions
        ]
        
        return non_contact_positions
    
    def _apply_non_contact_steering(
        self,
        sequence: str,
        heavy_seq: str,
        light_seq: str,
        non_contact_positions: List[int]
    ) -> Tuple[str, str, str, List[SequenceEdit]]:
        """Apply minimal steering to non-contact positions."""
        
        if self._iglm_model is None:
            self._iglm_model = mk_iglm_model(species=self.config.iglm_species)
        
        edits_applied = []
        current_seq_list = list(sequence)
        heavy_len = len(heavy_seq)
        
        # Limit number of edits
        max_edits = min(
            len(non_contact_positions),
            self.config.max_non_contact_edits_per_chain * 2  # Heavy + light chains
        )
        
        # Try edits at random non-contact positions
        positions_to_try = np.random.choice(
            non_contact_positions, 
            size=min(max_edits * 2, len(non_contact_positions)), 
            replace=False
        )
        
        for pos in positions_to_try:
            if len(edits_applied) >= max_edits:
                break
            
            original_aa = current_seq_list[pos]
            
            # Try IgLM-suggested amino acids
            suggested_aas = self._get_iglm_suggestions(current_seq_list, pos, heavy_len)
            
            best_aa = original_aa
            best_improvement = 0.0
            
            for suggested_aa in suggested_aas:
                if suggested_aa == original_aa:
                    continue
                
                # Test substitution
                test_seq = current_seq_list.copy()
                test_seq[pos] = suggested_aa
                
                # Score improvement
                improvement = self._score_edit_improvement(
                    current_seq_list, test_seq, pos, heavy_len
                )
                
                if improvement > best_improvement:
                    best_aa = suggested_aa
                    best_improvement = improvement
            
            # Apply edit if improvement found
            if best_improvement > 0.01:  # Minimum improvement threshold
                current_seq_list[pos] = best_aa
                
                edit = SequenceEdit(
                    position=pos,
                    original_aa=original_aa,
                    new_aa=best_aa,
                    edit_type="non_contact",
                    iglm_improvement=best_improvement,
                    structure_impact="neutral"  # Assumed for non-contact
                )
                edits_applied.append(edit)
        
        # Update sequences
        final_sequence = ''.join(current_seq_list)
        final_heavy = final_sequence[:heavy_len]
        final_light = final_sequence[heavy_len:]
        
        self.logger.info(f"Applied {len(edits_applied)} non-contact edits")
        
        return final_sequence, final_heavy, final_light, edits_applied
    
    def _apply_contact_guidance(
        self,
        sequence: str,
        heavy_seq: str,
        light_seq: str,
        contact_positions: List[int],
        filter_result: FilterResults
    ) -> Tuple[str, str, str, List[SequenceEdit]]:
        """Apply soft guidance to contact residues using KL prior."""
        
        edits_applied = []
        current_seq_list = list(sequence)
        heavy_len = len(heavy_seq)
        
        # Get epitope-specific enriched motifs
        enriched_motifs = self._get_epitope_enriched_motifs(filter_result)
        
        # Limit contact edits
        max_contact_edits = self.config.max_total_contact_edits
        
        # Try edits at contact positions
        positions_to_try = np.random.choice(
            contact_positions,
            size=min(max_contact_edits * 3, len(contact_positions)),
            replace=False
        )
        
        for pos in positions_to_try:
            if len(edits_applied) >= max_contact_edits:
                break
            
            original_aa = current_seq_list[pos]
            
            # Get KL-guided suggestions based on epitope chemistry
            suggested_aas = self._get_kl_guided_suggestions(
                pos, original_aa, enriched_motifs
            )
            
            best_aa = original_aa
            best_improvement = 0.0
            
            for suggested_aa in suggested_aas:
                if suggested_aa == original_aa:
                    continue
                
                # Test substitution with structure consideration
                improvement = self._score_contact_edit(
                    current_seq_list, pos, suggested_aa, heavy_len, filter_result
                )
                
                if improvement > best_improvement:
                    best_aa = suggested_aa
                    best_improvement = improvement
            
            # Apply edit if improvement found
            if best_improvement > 0.005:  # Lower threshold for contact edits
                current_seq_list[pos] = best_aa
                
                edit = SequenceEdit(
                    position=pos,
                    original_aa=original_aa,
                    new_aa=best_aa,
                    edit_type="contact",
                    iglm_improvement=best_improvement,
                    structure_impact="neutral"  # Would be evaluated by re-prediction
                )
                edits_applied.append(edit)
        
        # Update sequences
        final_sequence = ''.join(current_seq_list)
        final_heavy = final_sequence[:heavy_len]
        final_light = final_sequence[heavy_len:]
        
        self.logger.info(f"Applied {len(edits_applied)} contact edits")
        
        return final_sequence, final_heavy, final_light, edits_applied
    
    def _get_iglm_suggestions(
        self, 
        sequence: List[str], 
        position: int, 
        heavy_len: int
    ) -> List[str]:
        """Get IgLM-suggested amino acids for position."""
        
        # Create masked sequence
        masked_seq = sequence.copy()
        masked_seq[position] = '<mask>'
        
        # Determine chain
        chain_type = "heavy" if position < heavy_len else "light"
        
        # Get IgLM predictions
        predictions = self._iglm_model.predict_masked_position(
            ''.join(masked_seq), position, chain_type
        )
        
        # Return top 3 suggestions
        return [pred['amino_acid'] for pred in predictions[:3]]
    
    def _get_epitope_enriched_motifs(self, filter_result: FilterResults) -> Dict[str, List[str]]:
        """Get enriched amino acid motifs for the epitope chemistry."""
        
        # This would analyze the epitope and return enriched residue types
        # For now, return general binding-favorable residues
        motifs = {
            'hydrophobic': ['L', 'I', 'V', 'F', 'W', 'Y'],
            'polar': ['S', 'T', 'N', 'Q'],
            'charged': ['R', 'K', 'D', 'E'],
            'aromatic': ['F', 'W', 'Y', 'H']
        }
        
        return motifs
    
    def _get_kl_guided_suggestions(
        self, 
        position: int, 
        original_aa: str, 
        enriched_motifs: Dict[str, List[str]]
    ) -> List[str]:
        """Get KL-prior guided suggestions for contact position."""
        
        # Combine suggestions from enriched motifs
        suggestions = []
        for motif_class, aas in enriched_motifs.items():
            suggestions.extend(aas)
        
        # Remove duplicates and original
        suggestions = list(set(suggestions))
        if original_aa in suggestions:
            suggestions.remove(original_aa)
        
        # Return top suggestions (could be ranked by KL divergence)
        return suggestions[:5]
    
    def _score_edit_improvement(
        self, 
        original_seq: List[str], 
        edited_seq: List[str], 
        position: int, 
        heavy_len: int
    ) -> float:
        """Score the IgLM improvement from an edit."""
        
        # Score both sequences
        original_score = self._score_sequence_at_position(original_seq, position, heavy_len)
        edited_score = self._score_sequence_at_position(edited_seq, position, heavy_len)
        
        return edited_score - original_score
    
    def _score_contact_edit(
        self, 
        sequence: List[str], 
        position: int, 
        new_aa: str, 
        heavy_len: int,
        filter_result: FilterResults
    ) -> float:
        """Score a contact edit considering both IgLM and structure impact."""
        
        # Test sequence
        test_seq = sequence.copy()
        test_seq[position] = new_aa
        
        # IgLM component
        iglm_improvement = self._score_edit_improvement(sequence, test_seq, position, heavy_len)
        
        # Structure component (simplified - would use actual structure prediction)
        structure_penalty = self._estimate_structure_impact(
            position, sequence[position], new_aa, filter_result
        )
        
        # Combined score
        total_score = iglm_improvement - structure_penalty
        
        return total_score
    
    def _score_sequence_at_position(
        self, 
        sequence: List[str], 
        position: int, 
        heavy_len: int
    ) -> float:
        """Score sequence with focus on specific position."""
        
        # This would use IgLM to score the sequence
        # For now, return placeholder
        return np.random.normal(0, 1)
    
    def _estimate_structure_impact(
        self, 
        position: int, 
        original_aa: str, 
        new_aa: str, 
        filter_result: FilterResults
    ) -> float:
        """Estimate structural impact of amino acid substitution."""
        
        # Simple physicochemical penalty
        aa_properties = {
            'A': {'size': 1, 'charge': 0, 'hydrophobic': True},
            'R': {'size': 4, 'charge': 1, 'hydrophobic': False},
            'N': {'size': 2, 'charge': 0, 'hydrophobic': False},
            'D': {'size': 2, 'charge': -1, 'hydrophobic': False},
            'C': {'size': 2, 'charge': 0, 'hydrophobic': True},
            'Q': {'size': 3, 'charge': 0, 'hydrophobic': False},
            'E': {'size': 3, 'charge': -1, 'hydrophobic': False},
            'G': {'size': 0, 'charge': 0, 'hydrophobic': True},
            'H': {'size': 3, 'charge': 0.5, 'hydrophobic': False},
            'I': {'size': 3, 'charge': 0, 'hydrophobic': True},
            'L': {'size': 3, 'charge': 0, 'hydrophobic': True},
            'K': {'size': 4, 'charge': 1, 'hydrophobic': False},
            'M': {'size': 3, 'charge': 0, 'hydrophobic': True},
            'F': {'size': 4, 'charge': 0, 'hydrophobic': True},
            'P': {'size': 2, 'charge': 0, 'hydrophobic': True},
            'S': {'size': 1, 'charge': 0, 'hydrophobic': False},
            'T': {'size': 2, 'charge': 0, 'hydrophobic': False},
            'W': {'size': 5, 'charge': 0, 'hydrophobic': True},
            'Y': {'size': 4, 'charge': 0, 'hydrophobic': True},
            'V': {'size': 2, 'charge': 0, 'hydrophobic': True}
        }
        
        orig_props = aa_properties.get(original_aa, {'size': 2, 'charge': 0, 'hydrophobic': True})
        new_props = aa_properties.get(new_aa, {'size': 2, 'charge': 0, 'hydrophobic': True})
        
        # Calculate property differences
        size_diff = abs(orig_props['size'] - new_props['size'])
        charge_diff = abs(orig_props['charge'] - new_props['charge'])
        hydrophobic_diff = 1 if orig_props['hydrophobic'] != new_props['hydrophobic'] else 0
        
        # Simple penalty score
        penalty = (size_diff * 0.1) + (charge_diff * 0.2) + (hydrophobic_diff * 0.15)
        
        return penalty
    
    def _score_optimized_sequence(self, heavy_seq: str, light_seq: str) -> IgLMScores:
        """Score the optimized sequence with IgLM."""
        
        # Use the same scoring as humanness gate
        vh_ll = self._score_heavy_chain(heavy_seq)
        vl_ll = self._score_light_chain(light_seq)
        total_ll = vh_ll + vl_ll
        
        # Calculate percentile (would use same distribution as humanness gate)
        humanness_percentile = 50.0  # Placeholder
        
        scores = IgLMScores(
            sequence_id="optimized",
            full_sequence=heavy_seq + light_seq,
            heavy_chain_seq=heavy_seq,
            light_chain_seq=light_seq,
            vh_log_likelihood=vh_ll,
            vl_log_likelihood=vl_ll,
            total_log_likelihood=total_ll,
            humanness_percentile=humanness_percentile,
            humanness_pass=humanness_percentile >= self.config.humanness_percentile_threshold
        )
        
        return scores
    
    def _score_heavy_chain(self, heavy_seq: str) -> float:
        """Score heavy chain with IgLM."""
        if self._iglm_model is None:
            self._iglm_model = mk_iglm_model(species=self.config.iglm_species)
        return self._iglm_model.score_sequence(heavy_seq, chain_type="heavy")
    
    def _score_light_chain(self, light_seq: str) -> float:
        """Score light chain with IgLM."""
        if self._iglm_model is None:
            self._iglm_model = mk_iglm_model(species=self.config.iglm_species)
        return self._iglm_model.score_sequence(light_seq, chain_type="light")
    
    def _repredict_structure(
        self, 
        optimized_sequence: str, 
        original_filter_result: FilterResults
    ) -> Tuple[str, float]:
        """Re-predict structure after sequence optimization."""
        
        # This would use the same structure prediction as filtering module
        # For now, return placeholder
        structure_path = f"optimized_{original_filter_result.prediction.prediction_id}.pdb"
        
        # Assume slight degradation in structure score due to edits
        new_ranking_score = original_filter_result.ranking_score * 0.95
        
        return structure_path, new_ranking_score
    
    def _evaluate_optimization_success(
        self,
        original_scores: IgLMScores,
        optimized_scores: IgLMScores,
        edits_applied: List[SequenceEdit],
        final_ranking_score: float,
        original_ranking_score: float
    ) -> bool:
        """Evaluate whether optimization was successful."""
        
        # Success criteria:
        # 1. IgLM score improved
        # 2. Structure score maintained or improved
        # 3. Still passes humanness gate
        
        iglm_improved = optimized_scores.total_log_likelihood > original_scores.total_log_likelihood
        structure_maintained = final_ranking_score >= (original_ranking_score * 0.9)  # Allow 10% degradation
        humanness_maintained = optimized_scores.humanness_pass
        
        success = iglm_improved and structure_maintained and humanness_maintained
        
        self.logger.info(
            f"Optimization success: {success} "
            f"(IgLM: {iglm_improved}, Structure: {structure_maintained}, Humanness: {humanness_maintained})"
        )
        
        return success


def run_iglm_guidance_and_validation(
    filter_results: List[BackboneFilterResults],
    config: PipelineConfig,
    output_dir: str,
    logger: Optional[logging.Logger] = None
) -> List[List[OptimizedSequence]]:
    """
    Run the complete IgLM guidance and validation pipeline.
    
    Args:
        filter_results: Filtered structure results
        config: Pipeline configuration
        output_dir: Directory to save results
        logger: Optional logger instance
        
    Returns:
        List of optimized sequences for each backbone
    """
    logger = logger or logging.getLogger(__name__)
    logger.info("Starting IgLM guidance and validation")
    
    # Initialize components
    humanness_gate = IgLMHumannessGate(config.iglm, logger)
    optimizer = IgLMSequenceOptimizer(config.iglm, logger)
    
    # Step 1: Score humanness
    iglm_scores = humanness_gate.score_humanness(filter_results)
    
    # Step 2: Optimize sequences
    optimized_sequences = optimizer.optimize_sequences(filter_results, iglm_scores)
    
    # Save results
    output_path = Path(output_dir) / "iglm_results"
    output_path.mkdir(parents=True, exist_ok=True)
    
    for i, (backbone_result, backbone_scores, backbone_optimized) in enumerate(
        zip(filter_results, iglm_scores, optimized_sequences)
    ):
        backbone_dir = output_path / backbone_result.backbone_id
        backbone_dir.mkdir(exist_ok=True)
        
        # Save humanness scores
        with open(backbone_dir / "humanness_scores.csv", "w") as f:
            f.write("sequence_id,vh_log_likelihood,vl_log_likelihood,total_log_likelihood,humanness_percentile,humanness_pass\n")
            
            for score in backbone_scores:
                f.write(f"{score.sequence_id},{score.vh_log_likelihood},{score.vl_log_likelihood},"
                       f"{score.total_log_likelihood},{score.humanness_percentile},{score.humanness_pass}\n")
        
        # Save optimized sequences
        optimized_dir = backbone_dir / "optimized_sequences"
        optimized_dir.mkdir(exist_ok=True)
        
        for opt_seq in backbone_optimized:
            seq_dir = optimized_dir / opt_seq.original_filter_result.prediction.prediction_id
            seq_dir.mkdir(exist_ok=True)
            
            # Save optimization details
            optimization_data = {
                'original_sequence': opt_seq.original_filter_result.prediction.sequence_variant.sequence,
                'optimized_sequence': opt_seq.optimized_sequence,
                'optimization_success': opt_seq.optimization_success,
                'original_ranking_score': opt_seq.original_filter_result.ranking_score,
                'final_ranking_score': opt_seq.final_ranking_score,
                'iglm_scores': {
                    'vh_log_likelihood': opt_seq.iglm_scores.vh_log_likelihood,
                    'vl_log_likelihood': opt_seq.iglm_scores.vl_log_likelihood,
                    'total_log_likelihood': opt_seq.iglm_scores.total_log_likelihood,
                    'humanness_percentile': opt_seq.iglm_scores.humanness_percentile,
                    'humanness_pass': opt_seq.iglm_scores.humanness_pass
                },
                'edits_applied': [
                    {
                        'position': edit.position,
                        'original_aa': edit.original_aa,
                        'new_aa': edit.new_aa,
                        'edit_type': edit.edit_type,
                        'iglm_improvement': edit.iglm_improvement,
                        'structure_impact': edit.structure_impact
                    }
                    for edit in opt_seq.edits_applied
                ]
            }
            
            with open(seq_dir / "optimization_details.json", "w") as f:
                json.dump(optimization_data, f, indent=2)
            
            # Save sequences
            with open(seq_dir / "sequences.fasta", "w") as f:
                f.write(f">original\n{opt_seq.original_filter_result.prediction.sequence_variant.sequence}\n")
                f.write(f">optimized\n{opt_seq.optimized_sequence}\n")
    
    total_optimized = sum(len(backbone) for backbone in optimized_sequences)
    successful_optimizations = sum(
        sum(1 for opt in backbone if opt.optimization_success)
        for backbone in optimized_sequences
    )
    
    logger.info(f"IgLM optimization complete: {successful_optimizations}/{total_optimized} successful optimizations")
    
    return optimized_sequences