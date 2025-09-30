"""
CDR identification and labeling module

Based on: diffab-main/diffab/datasets/sabdab.py
"""

from typing import Optional, Dict, List
import torch
from .constants import CDR, ChothiaCDRRange


def label_heavy_chain_cdr(
    resseq: torch.Tensor,
    seq_map: Dict[str, int]
) -> torch.Tensor:
    """
    Label CDR regions in heavy chain using Chothia numbering.

    Args:
        resseq: Residue sequence numbers (Chothia numbering)
        seq_map: Mapping from sequential index to original index

    Returns:
        CDR labels tensor (0 for framework, CDR enum value for CDR regions)

    Based on: diffab/datasets/sabdab.py lines 74-104
    """
    cdr_flag = torch.zeros_like(resseq)

    for i_group, position in seq_map.items():
        if i_group not in ('H', 'A'):  # Heavy chain or antigen only
            continue

        resseq_this = resseq[position]
        cdr_flag_this = torch.zeros_like(resseq_this)

        for i_res, resseq_val in enumerate(resseq_this):
            resseq_val = resseq_val.item()

            # Check each CDR region
            if ChothiaCDRRange.H1[0] <= resseq_val <= ChothiaCDRRange.H1[1]:
                cdr_flag_this[i_res] = CDR.H1
            elif ChothiaCDRRange.H2[0] <= resseq_val <= ChothiaCDRRange.H2[1]:
                cdr_flag_this[i_res] = CDR.H2
            elif ChothiaCDRRange.H3[0] <= resseq_val <= ChothiaCDRRange.H3[1]:
                cdr_flag_this[i_res] = CDR.H3

        cdr_flag[position] = cdr_flag_this

    return cdr_flag


def label_light_chain_cdr(
    resseq: torch.Tensor,
    seq_map: Dict[str, int]
) -> torch.Tensor:
    """
    Label CDR regions in light chain using Chothia numbering.

    Args:
        resseq: Residue sequence numbers (Chothia numbering)
        seq_map: Mapping from sequential index to original index

    Returns:
        CDR labels tensor (0 for framework, CDR enum value for CDR regions)

    Based on: diffab/datasets/sabdab.py lines 107-134
    """
    cdr_flag = torch.zeros_like(resseq)

    for i_group, position in seq_map.items():
        if i_group not in ('L', 'A'):  # Light chain or antigen only
            continue

        resseq_this = resseq[position]
        cdr_flag_this = torch.zeros_like(resseq_this)

        for i_res, resseq_val in enumerate(resseq_this):
            resseq_val = resseq_val.item()

            # Check each CDR region
            if ChothiaCDRRange.L1[0] <= resseq_val <= ChothiaCDRRange.L1[1]:
                cdr_flag_this[i_res] = CDR.L1
            elif ChothiaCDRRange.L2[0] <= resseq_val <= ChothiaCDRRange.L2[1]:
                cdr_flag_this[i_res] = CDR.L2
            elif ChothiaCDRRange.L3[0] <= resseq_val <= ChothiaCDRRange.L3[1]:
                cdr_flag_this[i_res] = CDR.L3

        cdr_flag[position] = cdr_flag_this

    return cdr_flag


def extract_cdr_sequences(
    sequence: str,
    resseq: List[int],
    chain_type: str
) -> Dict[str, str]:
    """
    Extract CDR sequences from antibody chain.

    Args:
        sequence: Full amino acid sequence
        resseq: Residue sequence numbers (Chothia numbering)
        chain_type: 'H' for heavy, 'L' for light

    Returns:
        Dictionary mapping CDR names to sequences

    Example:
        >>> extract_cdr_sequences(seq, resseq, 'H')
        {'H1': 'GFTFSSYA', 'H2': 'ISGGS', 'H3': 'ARRGYLDY'}
    """
    assert chain_type in ('H', 'L'), f"Invalid chain type: {chain_type}"

    cdr_sequences = {}

    if chain_type == 'H':
        ranges = {
            'H1': ChothiaCDRRange.H1,
            'H2': ChothiaCDRRange.H2,
            'H3': ChothiaCDRRange.H3,
        }
    else:  # Light chain
        ranges = {
            'L1': ChothiaCDRRange.L1,
            'L2': ChothiaCDRRange.L2,
            'L3': ChothiaCDRRange.L3,
        }

    for cdr_name, (start, end) in ranges.items():
        cdr_residues = []
        for i, res_num in enumerate(resseq):
            if start <= res_num <= end:
                cdr_residues.append(sequence[i])

        if cdr_residues:
            cdr_sequences[cdr_name] = ''.join(cdr_residues)
        else:
            cdr_sequences[cdr_name] = ''

    return cdr_sequences


def validate_cdr_lengths(cdr_sequences: Dict[str, str]) -> Dict[str, bool]:
    """
    Validate CDR sequence lengths against typical ranges.

    Args:
        cdr_sequences: Dictionary mapping CDR names to sequences

    Returns:
        Dictionary mapping CDR names to validity flags

    Typical CDR length ranges:
        H1: 5-10 residues
        H2: 4-8 residues
        H3: 4-25 residues (most variable)
        L1: 6-15 residues
        L2: 3-7 residues
        L3: 4-12 residues
    """
    length_ranges = {
        'H1': (5, 10),
        'H2': (4, 8),
        'H3': (4, 25),
        'L1': (6, 15),
        'L2': (3, 7),
        'L3': (4, 12),
    }

    validation = {}
    for cdr_name, sequence in cdr_sequences.items():
        if cdr_name not in length_ranges:
            validation[cdr_name] = False
            continue

        min_len, max_len = length_ranges[cdr_name]
        length = len(sequence)
        validation[cdr_name] = min_len <= length <= max_len

    return validation


def get_cdr_info(
    heavy_seq: str,
    heavy_resseq: List[int],
    light_seq: str,
    light_resseq: List[int]
) -> Dict[str, any]:
    """
    Extract complete CDR information from antibody sequences.

    Args:
        heavy_seq: Heavy chain sequence
        heavy_resseq: Heavy chain residue numbers (Chothia)
        light_seq: Light chain sequence
        light_resseq: Light chain residue numbers (Chothia)

    Returns:
        Dictionary with CDR sequences and validation results
    """
    # Extract sequences
    heavy_cdrs = extract_cdr_sequences(heavy_seq, heavy_resseq, 'H')
    light_cdrs = extract_cdr_sequences(light_seq, light_resseq, 'L')

    # Combine all CDRs
    all_cdrs = {**heavy_cdrs, **light_cdrs}

    # Validate lengths
    validation = validate_cdr_lengths(all_cdrs)

    return {
        'sequences': all_cdrs,
        'validation': validation,
        'all_valid': all(validation.values()),
        'heavy_cdrs': heavy_cdrs,
        'light_cdrs': light_cdrs,
    }


def identify_paratope_residues(
    resseq: List[int],
    chain_type: str
) -> List[int]:
    """
    Identify potential paratope (antigen-binding) residues.

    The paratope typically includes all CDR residues plus some framework
    residues adjacent to CDRs.

    Args:
        resseq: Residue sequence numbers (Chothia numbering)
        chain_type: 'H' for heavy, 'L' for light

    Returns:
        List of indices corresponding to paratope residues
    """
    assert chain_type in ('H', 'L'), f"Invalid chain type: {chain_type}"

    paratope_indices = []

    if chain_type == 'H':
        ranges = [ChothiaCDRRange.H1, ChothiaCDRRange.H2, ChothiaCDRRange.H3]
    else:
        ranges = [ChothiaCDRRange.L1, ChothiaCDRRange.L2, ChothiaCDRRange.L3]

    for i, res_num in enumerate(resseq):
        for start, end in ranges:
            # Include CDR residues plus 2 residues on each side
            if (start - 2) <= res_num <= (end + 2):
                paratope_indices.append(i)
                break

    return paratope_indices