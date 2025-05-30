#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Phytovenomics Setup Entry Point

This module serves as the entry point for the phytovenomics-setup command,
which provides setup functionality for the Phytovenomics package.
"""

import sys
from phytovenomics.cli.setup_cli import main

if __name__ == '__main__':
    sys.exit(main())