"""
Integrated Antibody Design Pipeline
====================================

Complete end-to-end antibody design system combining:
- DiffAb: Structure preparation, CDR identification, renumbering
- Germinal: IgLM optimization, filtering, quality assessment
- RFdiffusion: Structure generation and design

Author: Phytovenomics Team
Version: 1.0.0
"""

__version__ = "1.0.0"

from antibody_pipeline.pipeline import AntibodyDesignPipeline

__all__ = ["AntibodyDesignPipeline"]