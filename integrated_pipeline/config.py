"""
Configuration management for the integrated pipeline.

This module provides configuration classes and utilities for managing
pipeline parameters, thresholds, and settings.
"""

import os
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Union
from pathlib import Path
import yaml


@dataclass
class GeometryConfig:
    """Configuration for pose and loop geometry generation."""
    
    # Framework settings
    framework_type: str = "human_igg1"  # Fix human IgG1 frameworks
    
    # Pose generation
    use_reference_complex: bool = False
    reference_complex_path: Optional[str] = None
    use_rigid_body_docking: bool = True
    epitope_cluster_constraint: bool = True
    
    # Loop generation methods
    use_rfdiffusion_loops: bool = True
    use_kic_loop_closure: bool = False
    use_ngk_loop_closure: bool = False
    use_equivariant_inverse_folding: bool = True
    
    # Diversity parameters
    min_diverse_loops: int = 50
    max_diverse_loops: int = 200
    diversity_rmsd_threshold: float = 2.0
    
    # CDR specific settings
    cdr_lengths: List[int] = field(default_factory=lambda: [10, 8, 12])  # H1, H2, H3
    cdr_positions: Optional[List[int]] = None


@dataclass
class SequenceConfig:
    """Configuration for sequence assignment on fixed backbones."""
    
    # MPNN settings
    use_antibody_mpnn: bool = True
    mpnn_model_path: str = "weights_abmpnn"
    condition_on_antigen: bool = True
    
    # Contact constraints
    constrain_contacts_to_cdrs: bool = True
    forbid_framework_contacts: bool = True
    
    # Sequence variants
    min_sequence_variants: int = 8
    max_sequence_variants: int = 16
    
    # Sampling parameters
    temperature: float = 1.0
    top_k: Optional[int] = None
    
    
@dataclass
class FilterConfig:
    """Configuration for complex prediction and filtering."""
    
    # Structure prediction
    use_af_multimer: bool = True
    use_af3: bool = False
    af_params_dir: str = "params"
    
    # Filter thresholds
    interface_pae_threshold: float = 5.0  # Median below 5 Å
    cdr_plddt_threshold: float = 0.8
    buried_surface_area_threshold: float = 800.0  # Å²
    paratope_fraction_threshold: float = 0.8  # 80% contacts from CDRs
    
    # Selection parameters
    top_n_per_backbone: int = 10
    
    
@dataclass
class IgLMConfig:
    """Configuration for IgLM guidance and validation."""
    
    # Humanness gate
    use_humanness_gate: bool = True
    humanness_percentile_threshold: float = 25.0  # Above 25th percentile
    therapeutic_distribution_path: Optional[str] = None
    
    # Minimal steering
    use_minimal_steering: bool = True
    steer_non_contact_positions: bool = True
    max_non_contact_edits_per_chain: int = 6
    
    # Soft guidance on contacts
    use_soft_contact_guidance: bool = True
    use_kl_prior: bool = True
    max_contact_edits_per_cdr: int = 3
    max_total_contact_edits: int = 3
    
    # Model settings
    iglm_model_path: Optional[str] = None
    iglm_species: str = "[HUMAN]"
    iglm_temperature: float = 1.0
    
    # Re-prediction after edits
    repredict_after_edits: bool = True
    keep_only_improved_poses: bool = True


@dataclass 
class ObjectiveConfig:
    """Configuration for broad vs specific objectives."""
    
    # Objective type
    objective_type: str = "broad"  # "broad" or "specific"
    
    # Broad binder settings
    average_across_cluster: bool = True
    add_variance_penalties: bool = True
    variance_penalty_weight: float = 0.1
    
    # Specific binder settings
    use_contrastive_decoys: bool = False
    require_unique_contacts: bool = False
    target_non_conserved_residues: bool = False
    

@dataclass
class RankingConfig:
    """Configuration for re-ranking and final selection."""
    
    # Ranking components
    use_independent_structure_score: bool = True
    use_interface_quality_score: bool = True
    use_paratope_fraction_score: bool = True
    use_rotamer_entropy_score: bool = True
    use_iglm_score: bool = True
    use_developability_score: bool = True
    
    # Weights for ranking components
    structure_weight: float = 0.25
    interface_weight: float = 0.20
    paratope_weight: float = 0.15
    rotamer_weight: float = 0.15
    iglm_weight: float = 0.15
    developability_weight: float = 0.10
    
    # Selection parameters
    final_selection_count: int = 50


@dataclass
class PipelineConfig:
    """Main configuration class for the integrated pipeline."""
    
    # Sub-configurations
    geometry: GeometryConfig = field(default_factory=GeometryConfig)
    sequence: SequenceConfig = field(default_factory=SequenceConfig)
    filtering: FilterConfig = field(default_factory=FilterConfig)
    iglm: IgLMConfig = field(default_factory=IgLMConfig)
    objective: ObjectiveConfig = field(default_factory=ObjectiveConfig)
    ranking: RankingConfig = field(default_factory=RankingConfig)
    
    # Global settings
    experiment_name: str = "phytovenomics_design"
    output_dir: str = "results"
    random_seed: int = 42
    
    # Input settings
    target_pdb_path: str = ""
    target_chain: str = "A"
    epitope_residues: List[int] = field(default_factory=list)
    epitope_clusters: Optional[Dict[str, List[int]]] = None
    
    # Resource settings
    use_gpu: bool = True
    max_gpu_memory: str = "40GB"
    num_parallel_workers: int = 1
    
    # Debug settings
    debug_mode: bool = False
    save_intermediate_results: bool = True
    verbose: bool = True
    
    @classmethod
    def from_yaml(cls, yaml_path: Union[str, Path]) -> 'PipelineConfig':
        """Load configuration from YAML file."""
        with open(yaml_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        
        # Convert nested dictionaries to dataclass instances
        if 'geometry' in config_dict:
            config_dict['geometry'] = GeometryConfig(**config_dict['geometry'])
        if 'sequence' in config_dict:
            config_dict['sequence'] = SequenceConfig(**config_dict['sequence'])
        if 'filtering' in config_dict:
            config_dict['filtering'] = FilterConfig(**config_dict['filtering'])
        if 'iglm' in config_dict:
            config_dict['iglm'] = IgLMConfig(**config_dict['iglm'])
        if 'objective' in config_dict:
            config_dict['objective'] = ObjectiveConfig(**config_dict['objective'])
        if 'ranking' in config_dict:
            config_dict['ranking'] = RankingConfig(**config_dict['ranking'])
            
        return cls(**config_dict)
    
    def to_yaml(self, yaml_path: Union[str, Path]) -> None:
        """Save configuration to YAML file."""
        # Convert dataclass instances to dictionaries
        config_dict = {
            'geometry': self.geometry.__dict__,
            'sequence': self.sequence.__dict__,
            'filtering': self.filtering.__dict__,
            'iglm': self.iglm.__dict__,
            'objective': self.objective.__dict__,
            'ranking': self.ranking.__dict__,
        }
        
        # Add global settings
        for key, value in self.__dict__.items():
            if key not in config_dict:
                config_dict[key] = value
                
        with open(yaml_path, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False, indent=2)
    
    def validate(self) -> None:
        """Validate configuration parameters."""
        errors = []
        
        # Check required paths
        if not self.target_pdb_path:
            errors.append("target_pdb_path is required")
        elif not os.path.exists(self.target_pdb_path):
            errors.append(f"target_pdb_path does not exist: {self.target_pdb_path}")
            
        # Check epitope residues
        if not self.epitope_residues and not self.epitope_clusters:
            errors.append("Either epitope_residues or epitope_clusters must be specified")
            
        # Check thresholds
        if not 0 <= self.filtering.cdr_plddt_threshold <= 1:
            errors.append("cdr_plddt_threshold must be between 0 and 1")
            
        if not 0 <= self.iglm.humanness_percentile_threshold <= 100:
            errors.append("humanness_percentile_threshold must be between 0 and 100")
            
        # Check ranking weights sum to 1
        total_weight = (
            self.ranking.structure_weight +
            self.ranking.interface_weight +
            self.ranking.paratope_weight +
            self.ranking.rotamer_weight +
            self.ranking.iglm_weight +
            self.ranking.developability_weight
        )
        if abs(total_weight - 1.0) > 0.01:
            errors.append(f"Ranking weights must sum to 1.0, got {total_weight}")
            
        if errors:
            raise ValueError("Configuration validation failed:\n" + "\n".join(errors))


def create_default_config() -> PipelineConfig:
    """Create a default configuration for testing."""
    config = PipelineConfig()
    config.target_pdb_path = "examples/target.pdb"
    config.epitope_residues = [25, 26, 27, 30, 31, 35]
    return config


def create_phytotoxin_config(target_name: str) -> PipelineConfig:
    """Create a configuration tailored for phytotoxin targets."""
    config = PipelineConfig()
    config.experiment_name = f"phytotoxin_{target_name}"
    
    # Adjust for phytotoxin-specific parameters
    config.geometry.max_diverse_loops = 150
    config.sequence.max_sequence_variants = 12
    config.filtering.buried_surface_area_threshold = 600.0  # Smaller toxins
    
    # More stringent humanness requirements
    config.iglm.humanness_percentile_threshold = 35.0
    config.iglm.max_contact_edits_per_cdr = 2
    
    # Broad binder objective for cross-reactivity
    config.objective.objective_type = "broad"
    config.objective.add_variance_penalties = True
    
    return config