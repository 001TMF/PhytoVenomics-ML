"""
Integrated DiffAb-Germinal Pipeline for Phytovenomics Antibody Design

This package integrates DiffAb and Germinal capabilities to implement a comprehensive
antibody design pipeline following the specified workflow:

1. Pose and loop geometry generation
2. Sequence assignment on fixed backbones  
3. Complex prediction and hard filters
4. IgLM guidance and validation
5. Re-ranking and final selection

The pipeline is designed for both broad and specific binder generation
against phytotoxin epitopes.
"""

__version__ = "1.0.0"
__author__ = "Phytovenomics-ML Pipeline"

from .pipeline import IntegratedPipeline
from .config import PipelineConfig

__all__ = ['IntegratedPipeline', 'PipelineConfig']