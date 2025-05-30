#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Phytovenomics Run Entry Point

This module serves as the entry point for the phytovenomics-run command,
which provides functionality for running various tools and analyses
within the Phytovenomics package.
"""

import sys
from phytovenomics.cli.run_cli import main

if __name__ == '__main__':
    sys.exit(main())