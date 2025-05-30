# phytovenomics/cli/__init__.py
"""
Phytovenomics CLI Package

This package contains command-line interfaces for the Phytovenomics project.
"""

from phytovenomics.cli.setup_cli import main as setup_main
from phytovenomics.cli.run_cli import main as run_main
from phytovenomics.cli.analyze_cli import main as analyze_main

__all__ = ["setup_main", "run_main", "analyze_main"]