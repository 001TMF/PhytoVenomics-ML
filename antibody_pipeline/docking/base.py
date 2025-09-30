"""
Base docking engine interface

Based on: diffab-main/diffab/tools/dock/base.py
"""

import abc
from typing import List


FilePath = str


class DockingEngine(abc.ABC):
    """
    Abstract base class for protein-protein docking engines.

    Implements context manager protocol for resource management.
    """

    @abc.abstractmethod
    def __enter__(self):
        """Enter context manager"""
        pass

    @abc.abstractmethod
    def __exit__(self, typ, value, traceback):
        """Exit context manager"""
        pass

    @abc.abstractmethod
    def set_receptor(self, pdb_path: FilePath):
        """
        Set receptor (antigen) structure.

        Args:
            pdb_path: Path to receptor PDB file
        """
        pass

    @abc.abstractmethod
    def set_ligand(self, pdb_path: FilePath):
        """
        Set ligand (antibody) structure.

        Args:
            pdb_path: Path to ligand PDB file
        """
        pass

    @abc.abstractmethod
    def dock(self) -> List[FilePath]:
        """
        Perform docking and return docked structures.

        Returns:
            List of paths to docked PDB files
        """
        pass