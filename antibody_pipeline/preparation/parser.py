"""
Structure Parsing Module

Based on:
- diffab-main/diffab/utils/protein/parsers.py
- diffab-main/diffab/datasets/custom.py
"""

import torch
import logging
from Bio.PDB import Selection, PDBParser
from Bio.PDB.Residue import Residue
from Bio.PDB import PDBExceptions
from typing import Dict, List, Tuple, Optional
from easydict import EasyDict

from antibody_pipeline.preparation.constants import (
    AA, max_num_heavyatoms,
    restype_to_heavyatom_names,
    BBHeavyAtom
)


class ParsingException(Exception):
    """Exception raised during structure parsing."""
    pass


def _get_residue_heavyatom_info(res: Residue) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Extract heavy atom coordinates and masks from residue.

    Args:
        res: BioPython Residue object

    Returns:
        Tuple of (positions, mask)
        positions: (max_num_heavyatoms, 3) tensor
        mask: (max_num_heavyatoms,) boolean tensor
    """
    pos_heavyatom = torch.zeros([max_num_heavyatoms, 3], dtype=torch.float)
    mask_heavyatom = torch.zeros([max_num_heavyatoms, ], dtype=torch.bool)

    restype = AA(res.get_resname())
    for idx, atom_name in enumerate(restype_to_heavyatom_names[restype]):
        if atom_name == '':
            continue
        if atom_name in res:
            pos_heavyatom[idx] = torch.tensor(
                res[atom_name].get_coord().tolist(),
                dtype=pos_heavyatom.dtype
            )
            mask_heavyatom[idx] = True

    return pos_heavyatom, mask_heavyatom


def parse_biopython_structure(
    entity,
    unknown_threshold: float = 1.0,
    max_resseq: Optional[int] = None
) -> Tuple[EasyDict, Dict]:
    """
    Parse BioPython structure entity into tensors.

    Based on: diffab-main/diffab/utils/protein/parsers.py

    Args:
        entity: BioPython entity (Chain, Model, or list of Chains)
        unknown_threshold: Maximum fraction of unknown residues allowed
        max_resseq: Maximum residue sequence number to include

    Returns:
        Tuple of (data, seq_map)
        data: EasyDict with structure tensors
        seq_map: Dictionary mapping (chain_id, resseq, icode) -> index
    """
    chains = Selection.unfold_entities(entity, 'C')
    chains.sort(key=lambda c: c.get_id())

    data = EasyDict({
        'chain_id': [],
        'resseq': [], 'icode': [], 'res_nb': [],
        'aa': [],
        'pos_heavyatom': [], 'mask_heavyatom': [],
    })

    tensor_types = {
        'resseq': torch.LongTensor,
        'res_nb': torch.LongTensor,
        'aa': torch.LongTensor,
        'pos_heavyatom': torch.stack,
        'mask_heavyatom': torch.stack,
    }

    count_aa, count_unk = 0, 0

    for i, chain in enumerate(chains):
        seq_this = 0  # Renumbering residues
        residues = Selection.unfold_entities(chain, 'R')
        residues.sort(key=lambda res: (res.get_id()[1], res.get_id()[2]))

        for _, res in enumerate(residues):
            resseq_this = int(res.get_id()[1])
            if max_resseq is not None and resseq_this > max_resseq:
                continue

            resname = res.get_resname()
            if not AA.is_aa(resname):
                continue
            if not (res.has_id('CA') and res.has_id('C') and res.has_id('N')):
                continue

            restype = AA(resname)
            count_aa += 1
            if restype == AA.UNK:
                count_unk += 1
                continue

            # Chain info
            data.chain_id.append(chain.get_id())

            # Residue types
            data.aa.append(restype)

            # Heavy atoms
            pos_heavyatom, mask_heavyatom = _get_residue_heavyatom_info(res)
            data.pos_heavyatom.append(pos_heavyatom)
            data.mask_heavyatom.append(mask_heavyatom)

            # Sequential number
            resseq_this = int(res.get_id()[1])
            icode_this = res.get_id()[2]
            if seq_this == 0:
                seq_this = 1
            else:
                d_CA_CA = torch.linalg.norm(
                    data.pos_heavyatom[-2][BBHeavyAtom.CA] - data.pos_heavyatom[-1][BBHeavyAtom.CA],
                    ord=2
                ).item()
                if d_CA_CA <= 4.0:
                    seq_this += 1
                else:
                    d_resseq = resseq_this - data.resseq[-1]
                    seq_this += max(2, d_resseq)

            data.resseq.append(resseq_this)
            data.icode.append(icode_this)
            data.res_nb.append(seq_this)

    if len(data.aa) == 0:
        raise ParsingException('No parsed residues.')

    if (count_unk / count_aa) >= unknown_threshold:
        raise ParsingException(
            f'Too many unknown residues, threshold {unknown_threshold:.2f}.'
        )

    seq_map = {}
    for i, (chain_id, resseq, icode) in enumerate(zip(data.chain_id, data.resseq, data.icode)):
        seq_map[(chain_id, resseq, icode)] = i

    for key, convert_fn in tensor_types.items():
        data[key] = convert_fn(data[key])

    return data, seq_map


def preprocess_antibody_structure(task: dict) -> Optional[dict]:
    """
    Preprocess antibody-antigen complex structure.

    Based on: diffab-main/diffab/datasets/custom.py

    Args:
        task: Dictionary with keys:
            - id: Structure identifier
            - pdb_path: Path to PDB file
            - heavy_id: Heavy chain identifier (default: 'H')
            - light_id: Light chain identifier (default: 'L')

    Returns:
        Dictionary with parsed structure data or None if parsing failed
    """
    from antibody_pipeline.preparation.cdr_identifier import label_heavy_chain_cdr, label_light_chain_cdr

    pdb_path = task['pdb_path']
    H_id = task.get('heavy_id', 'H')
    L_id = task.get('light_id', 'L')

    parser = PDBParser(QUIET=True)
    model = parser.get_structure('structure', pdb_path)[0]

    all_chain_ids = [c.id for c in model]

    parsed = {
        'id': task['id'],
        'heavy': None,
        'heavy_seqmap': None,
        'light': None,
        'light_seqmap': None,
        'antigen': None,
        'antigen_seqmap': None,
    }

    try:
        if H_id in all_chain_ids:
            (
                parsed['heavy'],
                parsed['heavy_seqmap']
            ) = label_heavy_chain_cdr(*parse_biopython_structure(
                model[H_id],
                max_resseq=113  # Chothia, end of Heavy chain Fv
            ))

        if L_id in all_chain_ids:
            (
                parsed['light'],
                parsed['light_seqmap']
            ) = label_light_chain_cdr(*parse_biopython_structure(
                model[L_id],
                max_resseq=106  # Chothia, end of Light chain Fv
            ))

        if parsed['heavy'] is None and parsed['light'] is None:
            raise ValueError(
                f'Neither valid antibody H-chain or L-chain is found. '
                f'Please ensure that the chain id of heavy chain is "{H_id}" '
                f'and the id of the light chain is "{L_id}".'
            )

        # Parse antigen chains (all non-antibody chains)
        ag_chain_ids = [cid for cid in all_chain_ids if cid not in (H_id, L_id)]
        if len(ag_chain_ids) > 0:
            chains = [model[c] for c in ag_chain_ids]
            (
                parsed['antigen'],
                parsed['antigen_seqmap']
            ) = parse_biopython_structure(chains)

    except (
        PDBExceptions.PDBConstructionException,
        ParsingException,
        KeyError,
        ValueError,
    ) as e:
        logging.warning('[{}] {}: {}'.format(
            task['id'],
            e.__class__.__name__,
            str(e)
        ))
        return None

    return parsed


class StructureParser:
    """
    High-level interface for structure parsing.
    """

    def parse_antibody_structure(
        self,
        pdb_path: str,
        structure_id: str = 'antibody',
        heavy_id: str = 'H',
        light_id: str = 'L'
    ) -> Optional[dict]:
        """
        Parse antibody-antigen complex.

        Args:
            pdb_path: Path to PDB file (should be Chothia numbered)
            structure_id: Structure identifier
            heavy_id: Heavy chain ID
            light_id: Light chain ID

        Returns:
            Parsed structure dictionary or None if failed
        """
        task = {
            'id': structure_id,
            'pdb_path': pdb_path,
            'heavy_id': heavy_id,
            'light_id': light_id,
        }

        return preprocess_antibody_structure(task)

    def parse_antigen_structure(
        self,
        pdb_path: str,
        structure_id: str = 'antigen'
    ) -> Optional[Tuple[EasyDict, Dict]]:
        """
        Parse antigen-only structure.

        Args:
            pdb_path: Path to PDB file
            structure_id: Structure identifier

        Returns:
            Tuple of (data, seq_map) or None if failed
        """
        try:
            parser = PDBParser(QUIET=True)
            model = parser.get_structure(structure_id, pdb_path)[0]
            return parse_biopython_structure(model)
        except Exception as e:
            logging.error(f"Failed to parse antigen: {e}")
            return None