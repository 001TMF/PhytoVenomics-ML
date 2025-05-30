#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Phytovenomics Analyze Entry Point

This module serves as the entry point for the phytovenomics-analyze command,
which provides analysis functionality for the Phytovenomics package.
"""

import sys
from phytovenomics.cli.analyze_cli import main

if __name__ == '__main__':
    sys.exit(main())