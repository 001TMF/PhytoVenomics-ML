"""
Main entry point for antibody_pipeline module

Allows running as: python -m antibody_pipeline
"""

from .cli import main
import sys

if __name__ == '__main__':
    sys.exit(main())