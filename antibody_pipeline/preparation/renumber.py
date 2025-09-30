"""
Antibody Renumbering Module

Based on: diffab-main/diffab/tools/renumber/run.py

Uses abnumber library for Chothia numbering scheme.
"""

import abnumber
from Bio import PDB
from Bio.PDB import Model, Chain, Residue, Selection
from Bio.Data import SCOPData
from typing import List, Tuple, Optional
import logging


def biopython_chain_to_sequence(chain: Chain.Chain) -> Tuple[str, List[Residue.Residue]]:
    """
    Extract sequence from BioPython chain.

    Args:
        chain: BioPython Chain object

    Returns:
        Tuple of (sequence, residue_list)
    """
    residue_list = Selection.unfold_entities(chain, 'R')
    seq = ''.join([SCOPData.protein_letters_3to1.get(r.resname, 'X') for r in residue_list])
    return seq, residue_list


def assign_chothia_numbering(seq: str) -> Tuple[List[Optional[Tuple[int, str]]], abnumber.Chain]:
    """
    Assign Chothia numbering to antibody sequence.

    Args:
        seq: Antibody sequence string

    Returns:
        Tuple of (numbering_list, abchain)
        numbering_list: List of (resseq, icode) tuples or None
        abchain: abnumber Chain object with metadata
    """
    abchain = abnumber.Chain(seq, scheme='chothia')
    offset = seq.index(abchain.seq)

    if not (offset >= 0):
        raise ValueError(
            'The identified Fv sequence is not a subsequence of the original sequence.'
        )

    numbers = [None for _ in range(len(seq))]
    for i, (pos, aa) in enumerate(abchain):
        resseq = pos.number
        icode = pos.letter if pos.letter else ' '
        numbers[i + offset] = (resseq, icode)

    return numbers, abchain


def renumber_biopython_chain(
    chain_id: str,
    residue_list: List[Residue.Residue],
    numbers: List[Optional[Tuple[int, str]]]
) -> Chain.Chain:
    """
    Create renumbered BioPython chain.

    Args:
        chain_id: Chain identifier
        residue_list: List of residue objects
        numbers: List of (resseq, icode) numbering tuples

    Returns:
        Renumbered Chain object
    """
    chain = Chain.Chain(chain_id)
    for residue, number in zip(residue_list, numbers):
        if number is None:
            continue
        residue = residue.copy()
        new_id = (residue.id[0], number[0], number[1])
        residue.id = new_id
        chain.add(residue)
    return chain


def renumber_antibody(
    in_pdb: str,
    out_pdb: str,
    return_other_chains: bool = False
) -> Tuple[List[str], List[str], Optional[List[str]]]:
    """
    Renumber antibody structure to Chothia scheme.

    Based on: diffab-main/diffab/tools/renumber/run.py

    Args:
        in_pdb: Input PDB file path
        out_pdb: Output PDB file path (renumbered)
        return_other_chains: Whether to return non-antibody chain IDs

    Returns:
        Tuple of (heavy_chains, light_chains, other_chains)
        heavy_chains: List of heavy chain IDs
        light_chains: List of light chain IDs
        other_chains: List of other chain IDs (if return_other_chains=True)
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(None, in_pdb)
    model = structure[0]
    model_new = Model.Model(0)

    heavy_chains, light_chains, other_chains = [], [], []

    for chain in model:
        try:
            seq, reslist = biopython_chain_to_sequence(chain)
            numbers, abchain = assign_chothia_numbering(seq)
            chain_new = renumber_biopython_chain(chain.id, reslist, numbers)

            logging.info(f'Renumbered chain {chain_new.id} ({abchain.chain_type})')

            if abchain.chain_type == 'H':
                heavy_chains.append(chain_new.id)
            elif abchain.chain_type in ('K', 'L'):
                light_chains.append(chain_new.id)
        except abnumber.ChainParseError as e:
            logging.info(f'Chain {chain.id} does not contain valid Fv: {str(e)}')
            chain_new = chain.copy()
            other_chains.append(chain_new.id)

        model_new.add(chain_new)

    # Save renumbered structure
    pdb_io = PDB.PDBIO()
    pdb_io.set_structure(model_new)
    pdb_io.save(out_pdb)

    if return_other_chains:
        return heavy_chains, light_chains, other_chains
    else:
        return heavy_chains, light_chains


class AntibodyRenumberer:
    """
    High-level interface for antibody renumbering.
    """

    def __init__(self, scheme: str = 'chothia'):
        """
        Initialize renumberer.

        Args:
            scheme: Numbering scheme (only 'chothia' supported currently)
        """
        self.scheme = scheme
        if scheme != 'chothia':
            raise ValueError("Only Chothia scheme is currently supported")

    def renumber(
        self,
        input_pdb: str,
        output_pdb: str
    ) -> dict:
        """
        Renumber antibody and return metadata.

        Args:
            input_pdb: Input PDB file
            output_pdb: Output renumbered PDB file

        Returns:
            Dictionary with chain information
        """
        heavy_chains, light_chains, other_chains = renumber_antibody(
            input_pdb,
            output_pdb,
            return_other_chains=True
        )

        return {
            'heavy_chains': heavy_chains,
            'light_chains': light_chains,
            'other_chains': other_chains,
            'scheme': self.scheme,
            'input_pdb': input_pdb,
            'output_pdb': output_pdb,
        }