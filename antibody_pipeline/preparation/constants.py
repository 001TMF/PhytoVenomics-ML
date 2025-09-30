"""
Protein constants and amino acid definitions

Based on: diffab-main/diffab/utils/protein/constants.py
"""

import torch
import enum
from typing import Dict


class CDR(enum.IntEnum):
    """CDR region identifiers"""
    H1 = 1
    H2 = 2
    H3 = 3
    L1 = 4
    L2 = 5
    L3 = 6


class ChothiaCDRRange:
    """
    Chothia CDR definitions.

    Based on: diffab-main/diffab/utils/protein/constants.py
    """
    H1 = (26, 32)
    H2 = (52, 56)
    H3 = (95, 102)

    L1 = (24, 34)
    L2 = (50, 56)
    L3 = (89, 97)

    @classmethod
    def to_cdr(cls, chain_type: str, resseq: int) -> CDR:
        """
        Determine CDR region from chain type and residue number.

        Args:
            chain_type: 'H' for heavy, 'L' for light
            resseq: Residue sequence number (Chothia numbering)

        Returns:
            CDR enum value or None
        """
        assert chain_type in ('H', 'L')

        if chain_type == 'H':
            if cls.H1[0] <= resseq <= cls.H1[1]:
                return CDR.H1
            elif cls.H2[0] <= resseq <= cls.H2[1]:
                return CDR.H2
            elif cls.H3[0] <= resseq <= cls.H3[1]:
                return CDR.H3
        elif chain_type == 'L':
            if cls.L1[0] <= resseq <= cls.L1[1]:
                return CDR.L1
            elif cls.L2[0] <= resseq <= cls.L2[1]:
                return CDR.L2
            elif cls.L3[0] <= resseq <= cls.L3[1]:
                return CDR.L3

        return None


class Fragment(enum.IntEnum):
    """Fragment type identifiers"""
    Heavy = 1
    Light = 2
    Antigen = 3


# Amino acid definitions
class AA(enum.IntEnum):
    """Amino acid types"""
    ALA = 0
    CYS = 1
    ASP = 2
    GLU = 3
    PHE = 4
    GLY = 5
    HIS = 6
    ILE = 7
    LYS = 8
    LEU = 9
    MET = 10
    ASN = 11
    PRO = 12
    GLN = 13
    ARG = 14
    SER = 15
    THR = 16
    VAL = 17
    TRP = 18
    TYR = 19
    UNK = 20

    @classmethod
    def is_aa(cls, resname: str) -> bool:
        """Check if residue name is a valid amino acid"""
        try:
            cls(resname)
            return True
        except (ValueError, KeyError):
            # Check if it's a non-standard residue
            return resname in non_standard_residue_substitutions

    def __call__(self, resname: str):
        """Get AA enum from residue name"""
        three_to_one = {
            'ALA': 0, 'CYS': 1, 'ASP': 2, 'GLU': 3, 'PHE': 4,
            'GLY': 5, 'HIS': 6, 'ILE': 7, 'LYS': 8, 'LEU': 9,
            'MET': 10, 'ASN': 11, 'PRO': 12, 'GLN': 13, 'ARG': 14,
            'SER': 15, 'THR': 16, 'VAL': 17, 'TRP': 18, 'TYR': 19,
        }

        # Handle non-standard residues
        if resname in non_standard_residue_substitutions:
            resname = non_standard_residue_substitutions[resname]

        if resname in three_to_one:
            return AA(three_to_one[resname])
        else:
            return AA.UNK


# Non-standard residue substitutions (from OpenMM)
non_standard_residue_substitutions = {
    '2AS': 'ASP', '3AH': 'HIS', '5HP': 'GLU', 'ACL': 'ARG', 'AGM': 'ARG',
    'AIB': 'ALA', 'ALM': 'ALA', 'ALO': 'THR', 'ALY': 'LYS', 'ARM': 'ARG',
    'ASA': 'ASP', 'ASB': 'ASP', 'ASK': 'ASP', 'ASL': 'ASP', 'ASQ': 'ASP',
    'MSE': 'MET', 'CSO': 'CYS', 'CME': 'CYS', # Common ones
}

# Backbone heavy atoms
class BBHeavyAtom(enum.IntEnum):
    N = 0
    CA = 1
    C = 2
    O = 3


# Heavy atom definitions
max_num_heavyatoms = 14

restype_to_heavyatom_names = {
    AA.ALA: ['N', 'CA', 'C', 'O', 'CB', '', '', '', '', '', '', '', '', ''],
    AA.CYS: ['N', 'CA', 'C', 'O', 'CB', 'SG', '', '', '', '', '', '', '', ''],
    AA.ASP: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2', '', '', '', '', '', ''],
    AA.GLU: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2', '', '', '', '', ''],
    AA.PHE: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', '', '', ''],
    AA.GLY: ['N', 'CA', 'C', 'O', '', '', '', '', '', '', '', '', '', ''],
    AA.HIS: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2', '', '', '', ''],
    AA.ILE: ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', '', '', '', '', '', ''],
    AA.LYS: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ', '', '', '', '', ''],
    AA.LEU: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', '', '', '', '', '', ''],
    AA.MET: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE', '', '', '', '', '', ''],
    AA.ASN: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2', '', '', '', '', '', ''],
    AA.PRO: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', '', '', '', '', '', '', ''],
    AA.GLN: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2', '', '', '', '', ''],
    AA.ARG: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2', '', '', ''],
    AA.SER: ['N', 'CA', 'C', 'O', 'CB', 'OG', '', '', '', '', '', '', '', ''],
    AA.THR: ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', '', '', '', '', '', '', ''],
    AA.VAL: ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', '', '', '', '', '', '', ''],
    AA.TRP: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    AA.TYR: ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH', '', ''],
    AA.UNK: ['', '', '', '', '', '', '', '', '', '', '', '', '', ''],
}

# Single letter to index mapping
ressymb_to_resindex = {
    'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4,
    'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9,
    'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14,
    'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19,
    'X': 20,
}

resindex_to_ressymb = {v: k for k, v in ressymb_to_resindex.items()}