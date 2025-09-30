"""
Complex Prediction and Hard Filters Module

This module implements the third stage of the pipeline:
1. Predict complexes with AF-Multimer or AF3
2. Filter by interface PAE, pLDDT on CDRs, buried surface area, and paratope fraction
3. Keep the top 5-10 per backbone

Integrates Germinal's filtering framework with structure prediction capabilities.
"""

import os
import numpy as np
import torch
import logging
from typing import List, Dict, Tuple, Optional, Union, Any
from pathlib import Path
from dataclasses import dataclass, field
import json

# Germinal filtering imports
from germinal.filters import filter_utils, af3, chai
from germinal.filters.filter_utils import run_filters
from germinal.utils.utils import calculate_clash_score

# ColabDesign imports for AF-Multimer
from colabdesign.af import mk_af_model

from .config import FilterConfig, PipelineConfig
from .sequence import BackboneDesigns, SequenceVariant


@dataclass
class StructurePrediction:
    """Container for structure prediction results."""
    prediction_id: str
    sequence_variant: SequenceVariant
    predicted_structure: str  # Path to PDB file
    confidence_scores: Dict[str, float]
    interface_metrics: Dict[str, float]
    quality_metrics: Dict[str, float]
    prediction_method: str  # "af_multimer" or "af3"
    
    
@dataclass
class FilterResults:
    """Container for filtering results on a structure prediction."""
    prediction: StructurePrediction
    filter_scores: Dict[str, float]
    filter_passes: Dict[str, bool]
    overall_pass: bool
    ranking_score: float
    

@dataclass
class BackboneFilterResults:
    """Container for all filter results for one backbone."""
    backbone_id: str
    all_predictions: List[StructurePrediction]
    filter_results: List[FilterResults]
    top_candidates: List[FilterResults]
    backbone_summary: Dict[str, Any]


class ComplexPredictor:
    """Structure prediction using AF-Multimer or AF3."""
    
    def __init__(self, config: FilterConfig, logger: Optional[logging.Logger] = None):
        self.config = config
        self.logger = logger or logging.getLogger(__name__)
        
        # Initialize models lazily
        self._af_multimer_model = None
        self._af3_model = None
        self._chai_model = None
        
    def predict_complexes(
        self,
        backbone_designs: List[BackboneDesigns],
        target_pdb: str
    ) -> List[List[StructurePrediction]]:
        """
        Predict structures for all sequence variants.
        
        Args:
            backbone_designs: List of backbone designs with sequence variants
            target_pdb: Path to target antigen PDB
            
        Returns:
            List of predictions for each backbone design
        """
        all_predictions = []
        
        for backbone_design in backbone_designs:
            self.logger.info(f"Predicting structures for backbone {backbone_design.backbone_id}")
            
            backbone_predictions = []
            
            for variant in backbone_design.sequence_variants:
                try:
                    prediction = self._predict_single_complex(
                        variant, target_pdb, backbone_design
                    )
                    backbone_predictions.append(prediction)
                    
                except Exception as e:
                    self.logger.error(f"Structure prediction failed for variant {variant.variant_id}: {e}")
                    continue
            
            all_predictions.append(backbone_predictions)
            
        return all_predictions
    
    def _predict_single_complex(
        self,
        variant: SequenceVariant,
        target_pdb: str,
        backbone_design: BackboneDesigns
    ) -> StructurePrediction:
        """Predict structure for a single sequence variant."""
        
        # Choose prediction method
        if self.config.use_af3:
            prediction = self._predict_with_af3(variant, target_pdb, backbone_design)
        elif self.config.use_af_multimer:
            prediction = self._predict_with_af_multimer(variant, target_pdb, backbone_design)
        else:
            raise ValueError("No structure prediction method enabled")
        
        return prediction
    
    def _predict_with_af_multimer(
        self,
        variant: SequenceVariant,
        target_pdb: str,
        backbone_design: BackboneDesigns
    ) -> StructurePrediction:
        """Predict structure using AlphaFold-Multimer."""
        
        if self._af_multimer_model is None:
            self._af_multimer_model = mk_af_model(
                protocol="complex",
                data_dir=self.config.af_params_dir,
                use_multimer=True
            )
        
        # Prepare input sequences
        target_sequence = self._extract_target_sequence(target_pdb)
        antibody_sequence = variant.sequence
        
        # Configure and run prediction
        self._af_multimer_model.prep_inputs(
            sequence=f"{target_sequence}:{antibody_sequence}",
            chain_linker=25  # Standard linker length
        )
        
        # Run prediction
        self._af_multimer_model.predict(
            num_models=1,
            num_recycles=3,
            save_best=True
        )
        
        # Save structure
        output_path = f"{variant.variant_id}_af_multimer.pdb"
        self._af_multimer_model.save_pdb(output_path)
        
        # Extract confidence and quality metrics
        confidence_scores = self._extract_af_multimer_confidence(self._af_multimer_model)
        interface_metrics = self._calculate_interface_metrics(output_path, target_pdb)
        quality_metrics = self._calculate_quality_metrics(output_path)
        
        prediction = StructurePrediction(
            prediction_id=f"{variant.variant_id}_af_multimer",
            sequence_variant=variant,
            predicted_structure=output_path,
            confidence_scores=confidence_scores,
            interface_metrics=interface_metrics,
            quality_metrics=quality_metrics,
            prediction_method="af_multimer"
        )
        
        return prediction
    
    def _predict_with_af3(
        self,
        variant: SequenceVariant,
        target_pdb: str,
        backbone_design: BackboneDesigns
    ) -> StructurePrediction:
        """Predict structure using AlphaFold3."""
        
        # Use Germinal's AF3 integration
        af3_config = {
            'target_chain': 'A',
            'binder_chain': 'B',
            'target_sequence': self._extract_target_sequence(target_pdb),
            'binder_sequence': variant.sequence
        }
        
        # Run AF3 prediction via Germinal
        af3_result = af3.run_af3_prediction(af3_config)
        
        # Extract metrics
        confidence_scores = af3_result.get('confidence_scores', {})
        interface_metrics = self._calculate_interface_metrics(
            af3_result['structure_path'], target_pdb
        )
        quality_metrics = self._calculate_quality_metrics(af3_result['structure_path'])
        
        prediction = StructurePrediction(
            prediction_id=f"{variant.variant_id}_af3",
            sequence_variant=variant,
            predicted_structure=af3_result['structure_path'],
            confidence_scores=confidence_scores,
            interface_metrics=interface_metrics,
            quality_metrics=quality_metrics,
            prediction_method="af3"
        )
        
        return prediction
    
    def _extract_af_multimer_confidence(self, af_model) -> Dict[str, float]:
        """Extract confidence scores from AF-Multimer model."""
        aux = af_model.aux
        
        # Extract key confidence metrics
        confidence_scores = {
            'plddt_mean': float(np.mean(aux['plddt'])),
            'plddt_interface': float(np.mean(aux['plddt'][-50:])),  # Last 50 residues (antibody)
            'ptm': float(aux.get('ptm', 0.0)),
            'iptm': float(aux.get('iptm', 0.0)),
            'ranking_confidence': float(aux.get('ranking_confidence', 0.0))
        }
        
        return confidence_scores
    
    def _calculate_interface_metrics(self, structure_pdb: str, target_pdb: str) -> Dict[str, float]:
        """Calculate interface-specific metrics."""
        
        # Use Germinal's utility functions
        interface_metrics = {
            'interface_pae': self._calculate_interface_pae(structure_pdb),
            'buried_surface_area': self._calculate_buried_surface_area(structure_pdb),
            'paratope_fraction': self._calculate_paratope_fraction(structure_pdb),
            'interface_contacts': self._count_interface_contacts(structure_pdb),
            'shape_complementarity': self._calculate_shape_complementarity(structure_pdb)
        }
        
        return interface_metrics
    
    def _calculate_quality_metrics(self, structure_pdb: str) -> Dict[str, float]:
        """Calculate overall structure quality metrics."""
        
        quality_metrics = {
            'clash_score': calculate_clash_score(structure_pdb, 2.4),
            'ca_clash_score': calculate_clash_score(structure_pdb, 2.5, only_ca=True),
            'ramachandran_favored': self._calculate_ramachandran_score(structure_pdb),
            'rotamer_outliers': self._calculate_rotamer_outliers(structure_pdb),
            'backbone_quality': self._assess_backbone_quality(structure_pdb)
        }
        
        return quality_metrics
    
    def _extract_target_sequence(self, target_pdb: str) -> str:
        """Extract amino acid sequence from target PDB."""
        # This would parse the PDB file and extract the sequence
        # For now, return placeholder
        return "PLACEHOLDER_TARGET_SEQUENCE"
    
    def _calculate_interface_pae(self, structure_pdb: str) -> float:
        """Calculate interface PAE."""
        # This would analyze PAE between chains
        return 3.5  # Placeholder
    
    def _calculate_buried_surface_area(self, structure_pdb: str) -> float:
        """Calculate buried surface area at interface."""
        # This would use geometric analysis
        return 850.0  # Placeholder
    
    def _calculate_paratope_fraction(self, structure_pdb: str) -> float:
        """Calculate fraction of contacts from CDRs (paratope)."""
        # This would identify CDR contacts vs framework contacts
        return 0.85  # Placeholder
    
    def _count_interface_contacts(self, structure_pdb: str) -> int:
        """Count contacts at the interface."""
        # This would identify residue-residue contacts
        return 25  # Placeholder
    
    def _calculate_shape_complementarity(self, structure_pdb: str) -> float:
        """Calculate shape complementarity at interface."""
        # This would use geometric analysis
        return 0.7  # Placeholder
        
    def _calculate_ramachandran_score(self, structure_pdb: str) -> float:
        """Calculate Ramachandran plot score."""
        return 0.95  # Placeholder
        
    def _calculate_rotamer_outliers(self, structure_pdb: str) -> float:
        """Calculate percentage of rotamer outliers."""
        return 0.02  # Placeholder
        
    def _assess_backbone_quality(self, structure_pdb: str) -> float:
        """Assess overall backbone quality."""
        return 0.9  # Placeholder


class StructureFilter:
    """Hard filters for structure quality and interface properties."""
    
    def __init__(self, config: FilterConfig, logger: Optional[logging.Logger] = None):
        self.config = config
        self.logger = logger or logging.getLogger(__name__)
        
    def apply_filters(
        self,
        predictions: List[List[StructurePrediction]]
    ) -> List[BackboneFilterResults]:
        """
        Apply hard filters to all structure predictions.
        
        Args:
            predictions: Nested list of predictions per backbone
            
        Returns:
            List of filter results per backbone
        """
        all_results = []
        
        for i, backbone_predictions in enumerate(predictions):
            backbone_id = f"backbone_{i}"
            
            self.logger.info(f"Applying filters to {len(backbone_predictions)} predictions for {backbone_id}")
            
            # Apply filters to each prediction
            filter_results = []
            for prediction in backbone_predictions:
                result = self._apply_single_filter(prediction)
                filter_results.append(result)
            
            # Select top candidates
            top_candidates = self._select_top_candidates(filter_results)
            
            # Create summary
            summary = self._create_backbone_summary(filter_results, top_candidates)
            
            backbone_result = BackboneFilterResults(
                backbone_id=backbone_id,
                all_predictions=backbone_predictions,
                filter_results=filter_results,
                top_candidates=top_candidates,
                backbone_summary=summary
            )
            
            all_results.append(backbone_result)
        
        return all_results
    
    def _apply_single_filter(self, prediction: StructurePrediction) -> FilterResults:
        """Apply all filters to a single prediction."""
        
        filter_scores = {}
        filter_passes = {}
        
        # Filter 1: Interface PAE
        interface_pae = prediction.interface_metrics['interface_pae']
        filter_scores['interface_pae'] = interface_pae
        filter_passes['interface_pae'] = interface_pae <= self.config.interface_pae_threshold
        
        # Filter 2: CDR pLDDT
        cdr_plddt = prediction.confidence_scores['plddt_interface']
        filter_scores['cdr_plddt'] = cdr_plddt
        filter_passes['cdr_plddt'] = cdr_plddt >= self.config.cdr_plddt_threshold
        
        # Filter 3: Buried surface area
        bsa = prediction.interface_metrics['buried_surface_area']
        filter_scores['buried_surface_area'] = bsa
        filter_passes['buried_surface_area'] = bsa >= self.config.buried_surface_area_threshold
        
        # Filter 4: Paratope fraction
        paratope_frac = prediction.interface_metrics['paratope_fraction']
        filter_scores['paratope_fraction'] = paratope_frac
        filter_passes['paratope_fraction'] = paratope_frac >= self.config.paratope_fraction_threshold
        
        # Filter 5: Clash score
        clash_score = prediction.quality_metrics['ca_clash_score']
        filter_scores['clash_score'] = clash_score
        filter_passes['clash_score'] = clash_score == 0  # No CA clashes
        
        # Filter 6: Interface contacts
        n_contacts = prediction.interface_metrics['interface_contacts']
        filter_scores['interface_contacts'] = n_contacts
        filter_passes['interface_contacts'] = n_contacts >= 8  # Minimum contacts
        
        # Overall pass/fail
        overall_pass = all(filter_passes.values())
        
        # Calculate ranking score for passed structures
        ranking_score = self._calculate_ranking_score(prediction, filter_scores) if overall_pass else 0.0
        
        result = FilterResults(
            prediction=prediction,
            filter_scores=filter_scores,
            filter_passes=filter_passes,
            overall_pass=overall_pass,
            ranking_score=ranking_score
        )
        
        return result
    
    def _calculate_ranking_score(
        self, 
        prediction: StructurePrediction, 
        filter_scores: Dict[str, float]
    ) -> float:
        """Calculate composite ranking score for passed structures."""
        
        # Normalize scores to 0-1 range and combine
        normalized_scores = {
            'confidence': prediction.confidence_scores['plddt_interface'],
            'interface_quality': 1.0 - (filter_scores['interface_pae'] / 10.0),  # Lower PAE is better
            'surface_area': min(filter_scores['buried_surface_area'] / 1500.0, 1.0),  # Cap at 1500
            'paratope': filter_scores['paratope_fraction'],
            'contacts': min(filter_scores['interface_contacts'] / 30.0, 1.0)  # Cap at 30
        }
        
        # Weighted combination
        weights = {
            'confidence': 0.3,
            'interface_quality': 0.25,
            'surface_area': 0.2,
            'paratope': 0.15,
            'contacts': 0.1
        }
        
        ranking_score = sum(
            weights[key] * normalized_scores[key] 
            for key in weights.keys()
        )
        
        return ranking_score
    
    def _select_top_candidates(self, filter_results: List[FilterResults]) -> List[FilterResults]:
        """Select top candidates that pass filters."""
        
        # Filter to only passing structures
        passing_results = [r for r in filter_results if r.overall_pass]
        
        if not passing_results:
            self.logger.warning("No structures passed all filters")
            return []
        
        # Sort by ranking score
        passing_results.sort(key=lambda x: x.ranking_score, reverse=True)
        
        # Take top N per backbone
        top_n = min(len(passing_results), self.config.top_n_per_backbone)
        top_candidates = passing_results[:top_n]
        
        self.logger.info(f"Selected {len(top_candidates)} top candidates from {len(filter_results)} total")
        
        return top_candidates
    
    def _create_backbone_summary(
        self, 
        all_results: List[FilterResults], 
        top_candidates: List[FilterResults]
    ) -> Dict[str, Any]:
        """Create summary statistics for the backbone."""
        
        summary = {
            'total_predictions': len(all_results),
            'passing_predictions': len([r for r in all_results if r.overall_pass]),
            'top_candidates': len(top_candidates),
            'pass_rate': len([r for r in all_results if r.overall_pass]) / len(all_results) if all_results else 0.0,
            'filter_pass_rates': {},
            'average_scores': {},
            'best_ranking_score': max([r.ranking_score for r in all_results]) if all_results else 0.0
        }
        
        # Calculate filter-specific pass rates
        if all_results:
            for filter_name in all_results[0].filter_passes.keys():
                passes = sum(1 for r in all_results if r.filter_passes[filter_name])
                summary['filter_pass_rates'][filter_name] = passes / len(all_results)
        
        # Calculate average scores for passed structures
        passed_results = [r for r in all_results if r.overall_pass]
        if passed_results:
            for score_name in passed_results[0].filter_scores.keys():
                scores = [r.filter_scores[score_name] for r in passed_results]
                summary['average_scores'][score_name] = np.mean(scores)
        
        return summary


def run_complex_prediction_and_filtering(
    backbone_designs: List[BackboneDesigns],
    config: PipelineConfig,
    output_dir: str,
    logger: Optional[logging.Logger] = None
) -> List[BackboneFilterResults]:
    """
    Run the complete complex prediction and filtering pipeline.
    
    Args:
        backbone_designs: List of backbone designs with sequence variants
        config: Pipeline configuration
        output_dir: Directory to save results
        logger: Optional logger instance
        
    Returns:
        List of filter results per backbone
    """
    logger = logger or logging.getLogger(__name__)
    logger.info("Starting complex prediction and filtering")
    
    # Initialize predictor and filter
    predictor = ComplexPredictor(config.filtering, logger)
    structure_filter = StructureFilter(config.filtering, logger)
    
    # Predict structures
    all_predictions = predictor.predict_complexes(
        backbone_designs=backbone_designs,
        target_pdb=config.target_pdb_path
    )
    
    # Apply filters
    filter_results = structure_filter.apply_filters(all_predictions)
    
    # Save results
    output_path = Path(output_dir) / "filtering_results"
    output_path.mkdir(parents=True, exist_ok=True)
    
    for backbone_result in filter_results:
        backbone_dir = output_path / backbone_result.backbone_id
        backbone_dir.mkdir(exist_ok=True)
        
        # Save backbone summary
        with open(backbone_dir / "summary.json", "w") as f:
            json.dump(backbone_result.backbone_summary, f, indent=2)
        
        # Save filter results
        results_dir = backbone_dir / "filter_results"
        results_dir.mkdir(exist_ok=True)
        
        with open(results_dir / "all_results.csv", "w") as f:
            f.write("prediction_id,overall_pass,ranking_score,interface_pae,cdr_plddt,buried_surface_area,paratope_fraction,clash_score,interface_contacts\n")
            
            for result in backbone_result.filter_results:
                f.write(f"{result.prediction.prediction_id},{result.overall_pass},{result.ranking_score},"
                       f"{result.filter_scores['interface_pae']},{result.filter_scores['cdr_plddt']},"
                       f"{result.filter_scores['buried_surface_area']},{result.filter_scores['paratope_fraction']},"
                       f"{result.filter_scores['clash_score']},{result.filter_scores['interface_contacts']}\n")
        
        # Save top candidates
        top_dir = backbone_dir / "top_candidates"
        top_dir.mkdir(exist_ok=True)
        
        for i, candidate in enumerate(backbone_result.top_candidates):
            candidate_dir = top_dir / f"rank_{i+1}_{candidate.prediction.prediction_id}"
            candidate_dir.mkdir(exist_ok=True)
            
            # Copy structure file
            import shutil
            shutil.copy(
                candidate.prediction.predicted_structure,
                candidate_dir / "structure.pdb"
            )
            
            # Save candidate details
            candidate_data = {
                'prediction_id': candidate.prediction.prediction_id,
                'sequence': candidate.prediction.sequence_variant.sequence,
                'ranking_score': candidate.ranking_score,
                'filter_scores': candidate.filter_scores,
                'filter_passes': candidate.filter_passes,
                'confidence_scores': candidate.prediction.confidence_scores,
                'interface_metrics': candidate.prediction.interface_metrics,
                'quality_metrics': candidate.prediction.quality_metrics
            }
            
            with open(candidate_dir / "details.json", "w") as f:
                json.dump(candidate_data, f, indent=2)
    
    total_passing = sum(len(br.top_candidates) for br in filter_results)
    logger.info(f"Filtering complete: {total_passing} structures passed all filters")
    
    return filter_results