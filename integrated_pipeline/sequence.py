"""
Sequence Assignment Module

This module implements the second stage of the pipeline:
1. Use antibody-tuned MPNN conditioned on antigen context
2. Constrain contacts to CDRs and forbid framework contacts  
3. Produce 8-16 sequence variants per backbone
4. Fixed backbone sequence design

Integrates Germinal's AbMPNN capabilities with the pipeline workflow.
"""

import os
import numpy as np
import torch
import logging
from typing import List, Dict, Tuple, Optional, Union
from pathlib import Path
from dataclasses import dataclass

# Germinal MPNN imports
from germinal.filters.redesign import run_abmpnn_redesign_pipeline
from colabdesign.mpnn import mk_mpnn_model

from .config import SequenceConfig, PipelineConfig
from .geometry import PoseGeometry, LoopGeometry


@dataclass
class SequenceVariant:
    """Container for a single sequence variant on a fixed backbone."""
    variant_id: str
    sequence: str
    heavy_chain_seq: str
    light_chain_seq: str
    mpnn_score: float
    sequence_recovery: float
    contact_recovery: float
    backbone_geometry: LoopGeometry
    antigen_context: str
    design_confidence: float


@dataclass 
class BackboneDesigns:
    """Container for all sequence variants on a single backbone."""
    backbone_id: str
    original_geometry: LoopGeometry
    sequence_variants: List[SequenceVariant]
    contact_map: np.ndarray
    interface_residues: List[int]
    framework_mask: np.ndarray  # True for framework positions


class SequenceAssigner:
    """Main class for sequence assignment on fixed backbones."""
    
    def __init__(self, config: SequenceConfig, logger: Optional[logging.Logger] = None):
        self.config = config
        self.logger = logger or logging.getLogger(__name__)
        
        # Initialize MPNN model lazily
        self._mpnn_model = None
        self._contact_predictor = None
        
    def assign_sequences_to_backbones(
        self,
        poses: List[PoseGeometry],
        target_pdb: str
    ) -> List[BackboneDesigns]:
        """
        Assign sequences to all backbone geometries using antibody-tuned MPNN.
        
        Args:
            poses: List of pose geometries with loop conformations
            target_pdb: Path to target antigen PDB
            
        Returns:
            List of BackboneDesigns with sequence variants
        """
        all_backbone_designs = []
        
        for pose in poses:
            self.logger.info(f"Processing pose {pose.pose_id} with {len(pose.loop_geometries)} loops")
            
            for loop_geometry in pose.loop_geometries:
                self.logger.info(f"Assigning sequences to loop {loop_geometry.loop_id}")
                
                # Create backbone design container
                backbone_design = self._create_backbone_design(
                    pose, loop_geometry, target_pdb
                )
                
                # Generate sequence variants
                variants = self._generate_sequence_variants(
                    backbone_design, target_pdb
                )
                
                backbone_design.sequence_variants = variants
                all_backbone_designs.append(backbone_design)
        
        self.logger.info(f"Generated sequences for {len(all_backbone_designs)} backbones")
        return all_backbone_designs
    
    def _create_backbone_design(
        self,
        pose: PoseGeometry,
        loop_geometry: LoopGeometry,
        target_pdb: str
    ) -> BackboneDesigns:
        """Create BackboneDesigns container with contact analysis."""
        
        # Create complex structure
        complex_pdb = self._create_complex_structure(pose, loop_geometry, target_pdb)
        
        # Analyze contacts and interface
        contact_map = self._analyze_contacts(complex_pdb)
        interface_residues = self._identify_interface_residues(contact_map)
        framework_mask = self._create_framework_mask(loop_geometry)
        
        backbone_design = BackboneDesigns(
            backbone_id=f"{pose.pose_id}_{loop_geometry.loop_id}",
            original_geometry=loop_geometry,
            sequence_variants=[],
            contact_map=contact_map,
            interface_residues=interface_residues,
            framework_mask=framework_mask
        )
        
        return backbone_design
    
    def _generate_sequence_variants(
        self,
        backbone_design: BackboneDesigns,
        target_pdb: str
    ) -> List[SequenceVariant]:
        """Generate sequence variants using antibody-tuned MPNN."""
        variants = []
        
        # Create structure file for MPNN
        structure_pdb = self._prepare_structure_for_mpnn(backbone_design)
        
        # Configure MPNN parameters
        mpnn_config = self._create_mpnn_config(backbone_design)
        
        # Run AbMPNN redesign pipeline
        try:
            mpnn_results = self._run_abmpnn_pipeline(
                structure_pdb, mpnn_config, target_pdb
            )
            
            # Convert MPNN results to SequenceVariant objects
            for i, result in enumerate(mpnn_results):
                variant = self._create_sequence_variant(
                    backbone_design, result, i
                )
                variants.append(variant)
                
        except Exception as e:
            self.logger.error(f"MPNN pipeline failed for backbone {backbone_design.backbone_id}: {e}")
            return []
        
        # Filter and validate variants
        validated_variants = self._validate_sequence_variants(variants, backbone_design)
        
        self.logger.info(f"Generated {len(validated_variants)} valid variants for backbone {backbone_design.backbone_id}")
        return validated_variants
    
    def _run_abmpnn_pipeline(
        self,
        structure_pdb: str,
        mpnn_config: Dict,
        target_pdb: str
    ) -> List[Dict]:
        """Run Germinal's AbMPNN redesign pipeline."""
        
        # Use Germinal's AbMPNN function
        abmpnn_sequences, success = run_abmpnn_redesign_pipeline(
            trajectory_pdb_af=structure_pdb,
            target_chain=mpnn_config['target_chain'],
            binder_chain=mpnn_config['binder_chain'],
            run_settings=mpnn_config['run_settings'],
            atom_distance_cutoff=mpnn_config['distance_cutoff']
        )
        
        if not success:
            raise RuntimeError("AbMPNN redesign pipeline failed")
        
        return abmpnn_sequences
    
    def _create_mpnn_config(self, backbone_design: BackboneDesigns) -> Dict:
        """Create configuration for MPNN run."""
        
        # Identify CDR positions to allow design
        cdr_positions = self._get_cdr_positions(backbone_design.original_geometry)
        
        # Create position mask - allow design only in CDRs
        design_mask = np.zeros(len(backbone_design.original_geometry.sequence), dtype=bool)
        for cdr_name, (start, end) in cdr_positions.items():
            design_mask[start:end+1] = True
        
        # Forbid framework contacts if configured
        if self.config.forbid_framework_contacts:
            framework_contacts = np.where(
                backbone_design.framework_mask & 
                (backbone_design.contact_map.sum(axis=1) > 0)
            )[0]
            design_mask[framework_contacts] = False
        
        config = {
            'target_chain': 'A',
            'binder_chain': 'B',
            'distance_cutoff': 5.0,
            'run_settings': {
                'mpnn_model_path': self.config.mpnn_model_path,
                'temperature': self.config.temperature,
                'num_sequences': self.config.max_sequence_variants,
                'design_mask': design_mask.tolist(),
                'condition_on_antigen': self.config.condition_on_antigen,
                'constrain_contacts_to_cdrs': self.config.constrain_contacts_to_cdrs,
                'top_k': self.config.top_k
            }
        }
        
        return config
    
    def _create_sequence_variant(
        self,
        backbone_design: BackboneDesigns,
        mpnn_result: Dict,
        variant_index: int
    ) -> SequenceVariant:
        """Create SequenceVariant from MPNN result."""
        
        full_sequence = mpnn_result['seq']
        
        # Split into heavy and light chains based on geometry
        heavy_chain_len = len([pos for cdr, (start, end) in 
                              backbone_design.original_geometry.cdr_assignments.items() 
                              if cdr.startswith('H')])
        
        heavy_chain_seq = full_sequence[:heavy_chain_len]
        light_chain_seq = full_sequence[heavy_chain_len:]
        
        # Calculate metrics
        sequence_recovery = self._calculate_sequence_recovery(
            full_sequence, backbone_design.original_geometry.sequence
        )
        
        contact_recovery = self._calculate_contact_recovery(
            backbone_design, full_sequence
        )
        
        variant = SequenceVariant(
            variant_id=f"{backbone_design.backbone_id}_var_{variant_index}",
            sequence=full_sequence,
            heavy_chain_seq=heavy_chain_seq,
            light_chain_seq=light_chain_seq,
            mpnn_score=mpnn_result.get('score', 0.0),
            sequence_recovery=sequence_recovery,
            contact_recovery=contact_recovery,
            backbone_geometry=backbone_design.original_geometry,
            antigen_context=self._extract_antigen_context(backbone_design),
            design_confidence=mpnn_result.get('confidence', 0.0)
        )
        
        return variant
    
    def _validate_sequence_variants(
        self,
        variants: List[SequenceVariant],
        backbone_design: BackboneDesigns
    ) -> List[SequenceVariant]:
        """Validate and filter sequence variants."""
        valid_variants = []
        
        for variant in variants:
            # Check minimum quality thresholds
            if variant.mpnn_score < -10.0:  # Very poor MPNN score
                continue
                
            if variant.sequence_recovery < 0.1:  # Too different from original
                continue
                
            # Check for problematic sequences
            if self._has_problematic_motifs(variant.sequence):
                continue
                
            # Check CDR contact constraints
            if (self.config.constrain_contacts_to_cdrs and 
                not self._validates_cdr_contacts(variant, backbone_design)):
                continue
                
            valid_variants.append(variant)
        
        # Sort by MPNN score and limit number
        valid_variants.sort(key=lambda x: x.mpnn_score, reverse=True)
        
        # Ensure we have minimum and maximum variants
        if len(valid_variants) < self.config.min_sequence_variants:
            self.logger.warning(
                f"Only {len(valid_variants)} valid variants for backbone "
                f"{backbone_design.backbone_id}, minimum is {self.config.min_sequence_variants}"
            )
        
        return valid_variants[:self.config.max_sequence_variants]
    
    def _analyze_contacts(self, complex_pdb: str) -> np.ndarray:
        """Analyze contacts in the antibody-antigen complex."""
        # This would use structural analysis to identify contacts
        # For now, return placeholder contact map
        n_residues = 200  # Placeholder
        contact_map = np.random.rand(n_residues, n_residues) > 0.95
        return contact_map.astype(float)
    
    def _identify_interface_residues(self, contact_map: np.ndarray) -> List[int]:
        """Identify residues at the antibody-antigen interface."""
        # Residues with any inter-chain contacts
        interface_mask = contact_map.sum(axis=1) > 0
        return np.where(interface_mask)[0].tolist()
    
    def _create_framework_mask(self, loop_geometry: LoopGeometry) -> np.ndarray:
        """Create mask indicating framework vs CDR positions."""
        n_residues = len(loop_geometry.sequence)
        framework_mask = np.ones(n_residues, dtype=bool)
        
        # Mark CDR positions as False (not framework)
        for cdr_name, (start, end) in loop_geometry.cdr_assignments.items():
            framework_mask[start:end+1] = False
            
        return framework_mask
    
    def _get_cdr_positions(self, loop_geometry: LoopGeometry) -> Dict[str, Tuple[int, int]]:
        """Get CDR position ranges."""
        return loop_geometry.cdr_assignments
    
    def _calculate_sequence_recovery(self, new_seq: str, original_seq: str) -> float:
        """Calculate sequence recovery between new and original sequences."""
        if len(new_seq) != len(original_seq):
            return 0.0
        
        matches = sum(1 for a, b in zip(new_seq, original_seq) if a == b)
        return matches / len(new_seq)
    
    def _calculate_contact_recovery(
        self, 
        backbone_design: BackboneDesigns, 
        new_sequence: str
    ) -> float:
        """Calculate recovery of important contacts."""
        # This would analyze how well the new sequence maintains key contacts
        # For now, return placeholder value
        return 0.8
    
    def _extract_antigen_context(self, backbone_design: BackboneDesigns) -> str:
        """Extract relevant antigen context for the design."""
        # This would extract antigen sequence/structure context
        # For now, return placeholder
        return "ANTIGEN_CONTEXT"
    
    def _has_problematic_motifs(self, sequence: str) -> bool:
        """Check for problematic sequence motifs."""
        problematic_motifs = [
            "NG",    # Deamidation site
            "DG",    # Isomerization site  
            "CW",    # Tryptophan oxidation near cysteine
            "GGGG",  # Flexible linker (problematic in CDRs)
        ]
        
        return any(motif in sequence for motif in problematic_motifs)
    
    def _validates_cdr_contacts(
        self, 
        variant: SequenceVariant, 
        backbone_design: BackboneDesigns
    ) -> bool:
        """Validate that contacts are constrained to CDRs."""
        if not self.config.constrain_contacts_to_cdrs:
            return True
        
        # Check that framework residues don't make contacts
        framework_interface = (
            backbone_design.framework_mask & 
            np.isin(np.arange(len(backbone_design.framework_mask)), 
                   backbone_design.interface_residues)
        )
        
        # Allow some tolerance for framework contacts
        framework_contact_fraction = framework_interface.sum() / len(backbone_design.interface_residues)
        
        return framework_contact_fraction < 0.2  # Less than 20% framework contacts
    
    def _create_complex_structure(
        self, 
        pose: PoseGeometry, 
        loop_geometry: LoopGeometry, 
        target_pdb: str
    ) -> str:
        """Create complete antibody-antigen complex structure."""
        # This would combine antigen, framework, and loop coordinates
        # Return path to temporary PDB file
        return "complex_with_loops.pdb"
    
    def _prepare_structure_for_mpnn(self, backbone_design: BackboneDesigns) -> str:
        """Prepare structure file for MPNN input."""
        # This would format the structure appropriately for MPNN
        # Return path to MPNN-ready PDB file
        return "mpnn_input.pdb"


def run_sequence_assignment(
    poses: List[PoseGeometry],
    config: PipelineConfig,
    output_dir: str,
    logger: Optional[logging.Logger] = None
) -> List[BackboneDesigns]:
    """
    Run the complete sequence assignment pipeline.
    
    Args:
        poses: List of pose geometries from geometry generation
        config: Pipeline configuration
        output_dir: Directory to save results
        logger: Optional logger instance
        
    Returns:
        List of backbone designs with sequence variants
    """
    logger = logger or logging.getLogger(__name__)
    logger.info("Starting sequence assignment on fixed backbones")
    
    # Initialize assigner
    assigner = SequenceAssigner(config.sequence, logger)
    
    # Assign sequences to all backbones
    backbone_designs = assigner.assign_sequences_to_backbones(
        poses=poses,
        target_pdb=config.target_pdb_path
    )
    
    # Save results
    output_path = Path(output_dir) / "sequence_results"
    output_path.mkdir(parents=True, exist_ok=True)
    
    for backbone_design in backbone_designs:
        backbone_dir = output_path / backbone_design.backbone_id
        backbone_dir.mkdir(exist_ok=True)
        
        # Save contact analysis
        np.save(backbone_dir / "contact_map.npy", backbone_design.contact_map)
        np.save(backbone_dir / "framework_mask.npy", backbone_design.framework_mask)
        
        with open(backbone_dir / "interface_residues.txt", "w") as f:
            f.write(",".join(map(str, backbone_design.interface_residues)))
        
        # Save sequence variants
        variants_dir = backbone_dir / "variants"
        variants_dir.mkdir(exist_ok=True)
        
        with open(variants_dir / "summary.csv", "w") as f:
            f.write("variant_id,sequence,heavy_chain_seq,light_chain_seq,mpnn_score,sequence_recovery,contact_recovery,design_confidence\n")
            
            for variant in backbone_design.sequence_variants:
                f.write(f"{variant.variant_id},{variant.sequence},{variant.heavy_chain_seq},"
                       f"{variant.light_chain_seq},{variant.mpnn_score},{variant.sequence_recovery},"
                       f"{variant.contact_recovery},{variant.design_confidence}\n")
                
                # Save individual variant
                variant_file = variants_dir / f"{variant.variant_id}.txt"
                with open(variant_file, "w") as vf:
                    vf.write(f"Sequence: {variant.sequence}\n")
                    vf.write(f"Heavy Chain: {variant.heavy_chain_seq}\n")
                    vf.write(f"Light Chain: {variant.light_chain_seq}\n")
                    vf.write(f"MPNN Score: {variant.mpnn_score}\n")
                    vf.write(f"Sequence Recovery: {variant.sequence_recovery}\n")
                    vf.write(f"Contact Recovery: {variant.contact_recovery}\n")
                    vf.write(f"Antigen Context: {variant.antigen_context}\n")
    
    total_variants = sum(len(bd.sequence_variants) for bd in backbone_designs)
    logger.info(f"Generated {total_variants} sequence variants across {len(backbone_designs)} backbones")
    
    return backbone_designs