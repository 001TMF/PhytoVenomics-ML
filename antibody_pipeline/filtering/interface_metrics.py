"""
Antibody-antigen interface metrics

Based on: germinal-main/germinal/filters/filter_utils.py
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from Bio import PDB


class InterfaceMetrics:
    """
    Calculate interface metrics for antibody-antigen complexes.

    Metrics include:
    - Binding energy (dG)
    - Solvent accessible surface area (SASA)
    - Interface contacts
    - Hydrogen bonds
    - pDockQ (predicted DockQ score)
    """

    def __init__(self, contact_distance: float = 5.0):
        """
        Initialize interface metrics calculator.

        Args:
            contact_distance: Distance threshold for contacts (Angstroms)
        """
        self.contact_distance = contact_distance

    def calculate_contacts(
        self,
        pdb_path: str,
        chain1: str,
        chain2: str
    ) -> Tuple[int, List[Tuple[str, str]]]:
        """
        Calculate number of residue contacts between two chains.

        Args:
            pdb_path: Path to PDB file
            chain1: First chain ID
            chain2: Second chain ID

        Returns:
            Tuple of (num_contacts, contact_list)
        """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("complex", pdb_path)

        # Get residues from each chain
        residues1 = []
        residues2 = []

        for model in structure:
            if chain1 in model:
                residues1 = list(model[chain1].get_residues())
            if chain2 in model:
                residues2 = list(model[chain2].get_residues())

        # Calculate contacts
        contacts = []
        for res1 in residues1:
            for res2 in residues2:
                # Check if any atoms are within contact distance
                for atom1 in res1:
                    for atom2 in res2:
                        distance = atom1 - atom2
                        if distance <= self.contact_distance:
                            res1_id = f"{res1.get_resname()}{res1.id[1]}"
                            res2_id = f"{res2.get_resname()}{res2.id[1]}"
                            contacts.append((res1_id, res2_id))
                            break
                    else:
                        continue
                    break

        return len(contacts), contacts

    def calculate_interface_residues(
        self,
        pdb_path: str,
        antibody_chain: str,
        antigen_chain: str
    ) -> Dict[str, List[int]]:
        """
        Identify interface residues in antibody and antigen.

        Args:
            pdb_path: Path to PDB file
            antibody_chain: Antibody chain ID
            antigen_chain: Antigen chain ID

        Returns:
            Dictionary with antibody and antigen interface residues
        """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("complex", pdb_path)

        antibody_residues = []
        antigen_residues = []

        for model in structure:
            if antibody_chain in model and antigen_chain in model:
                ab_chain = model[antibody_chain]
                ag_chain = model[antigen_chain]

                # Find interface residues
                for ab_res in ab_chain:
                    for ag_res in ag_chain:
                        min_dist = float('inf')
                        for ab_atom in ab_res:
                            for ag_atom in ag_res:
                                dist = ab_atom - ag_atom
                                if dist < min_dist:
                                    min_dist = dist

                        if min_dist <= self.contact_distance:
                            if ab_res.id[1] not in antibody_residues:
                                antibody_residues.append(ab_res.id[1])
                            if ag_res.id[1] not in antigen_residues:
                                antigen_residues.append(ag_res.id[1])

        return {
            'antibody_interface': sorted(antibody_residues),
            'antigen_interface': sorted(antigen_residues),
            'num_antibody_interface': len(antibody_residues),
            'num_antigen_interface': len(antigen_residues),
        }

    def calculate_cdr_interface_percentage(
        self,
        pdb_path: str,
        antibody_chain: str,
        antigen_chain: str,
        cdr_positions: List[int]
    ) -> float:
        """
        Calculate percentage of interface residues that are CDR residues.

        Args:
            pdb_path: Path to PDB file
            antibody_chain: Antibody chain ID
            antigen_chain: Antigen chain ID
            cdr_positions: List of CDR residue positions

        Returns:
            Percentage of interface residues in CDRs
        """
        interface_residues = self.calculate_interface_residues(
            pdb_path, antibody_chain, antigen_chain
        )

        ab_interface = interface_residues['antibody_interface']

        # Count how many interface residues are in CDRs
        cdr_interface_count = sum(
            1 for res in ab_interface if res in cdr_positions
        )

        if len(ab_interface) == 0:
            return 0.0

        return 100.0 * cdr_interface_count / len(ab_interface)

    def calculate_hotspot_contacts(
        self,
        pdb_path: str,
        antibody_chain: str,
        antigen_chain: str,
        hotspot_residues: List[int],
        cdr3_positions: Optional[List[int]] = None
    ) -> Dict[str, any]:
        """
        Calculate contacts between antibody and antigen hotspot residues.

        Args:
            pdb_path: Path to PDB file
            antibody_chain: Antibody chain ID
            antigen_chain: Antigen chain ID
            hotspot_residues: List of hotspot residue positions in antigen
            cdr3_positions: Optional CDR3 residue positions

        Returns:
            Dictionary with hotspot contact metrics
        """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("complex", pdb_path)

        hotspot_contacts = []
        cdr3_hotspot_contacts = []

        for model in structure:
            if antibody_chain in model and antigen_chain in model:
                ab_chain = model[antibody_chain]
                ag_chain = model[antigen_chain]

                # Check contacts with hotspots
                for ab_res in ab_chain:
                    ab_pos = ab_res.id[1]
                    is_cdr3 = cdr3_positions and (ab_pos in cdr3_positions)

                    for ag_res in ag_chain:
                        ag_pos = ag_res.id[1]

                        if ag_pos in hotspot_residues:
                            # Check if in contact
                            min_dist = float('inf')
                            for ab_atom in ab_res:
                                for ag_atom in ag_res:
                                    dist = ab_atom - ag_atom
                                    if dist < min_dist:
                                        min_dist = dist

                            if min_dist <= self.contact_distance:
                                hotspot_contacts.append((ab_pos, ag_pos))
                                if is_cdr3:
                                    cdr3_hotspot_contacts.append((ab_pos, ag_pos))

        return {
            'num_hotspot_contacts': len(hotspot_contacts),
            'num_cdr3_hotspot_contacts': len(cdr3_hotspot_contacts),
            'hotspot_contacts': hotspot_contacts,
            'cdr3_hotspot_contacts': cdr3_hotspot_contacts,
        }

    def calculate_buried_surface_area(
        self,
        complex_pdb: str,
        antibody_pdb: str,
        antigen_pdb: str
    ) -> Dict[str, float]:
        """
        Calculate buried surface area upon binding.

        Args:
            complex_pdb: Path to complex PDB
            antibody_pdb: Path to isolated antibody PDB
            antigen_pdb: Path to isolated antigen PDB

        Returns:
            Dictionary with BSA metrics
        """
        try:
            from Bio.PDB import SASA

            parser = PDB.PDBParser(QUIET=True)

            # Calculate SASA for complex
            complex_struct = parser.get_structure("complex", complex_pdb)
            complex_sasa = SASA.ShrakeRupley()
            complex_sasa.compute(complex_struct, level="R")
            complex_area = sum(
                residue.sasa for model in complex_struct
                for chain in model for residue in chain
            )

            # Calculate SASA for antibody alone
            ab_struct = parser.get_structure("antibody", antibody_pdb)
            ab_sasa = SASA.ShrakeRupley()
            ab_sasa.compute(ab_struct, level="R")
            ab_area = sum(
                residue.sasa for model in ab_struct
                for chain in model for residue in chain
            )

            # Calculate SASA for antigen alone
            ag_struct = parser.get_structure("antigen", antigen_pdb)
            ag_sasa = SASA.ShrakeRupley()
            ag_sasa.compute(ag_struct, level="R")
            ag_area = sum(
                residue.sasa for model in ag_struct
                for chain in model for residue in chain
            )

            # BSA = (SASA_ab + SASA_ag - SASA_complex) / 2
            bsa = (ab_area + ag_area - complex_area) / 2.0

            return {
                'bsa': bsa,
                'antibody_sasa': ab_area,
                'antigen_sasa': ag_area,
                'complex_sasa': complex_area,
            }

        except Exception as e:
            print(f"Warning: Could not calculate BSA: {e}")
            return {
                'bsa': 0.0,
                'antibody_sasa': 0.0,
                'antigen_sasa': 0.0,
                'complex_sasa': 0.0,
            }

    def calculate_pdockq(
        self,
        pdb_path: str,
        pae_matrix: Optional[np.ndarray] = None,
        plddt_values: Optional[np.ndarray] = None
    ) -> float:
        """
        Calculate predicted DockQ score.

        pDockQ estimates the quality of protein-protein docking
        using predicted aligned error (PAE) and pLDDT.

        Args:
            pdb_path: Path to PDB file
            pae_matrix: Predicted aligned error matrix (optional)
            plddt_values: Per-residue pLDDT values (optional)

        Returns:
            pDockQ score (0-1, higher is better)
        """
        # Simplified pDockQ calculation
        # Full implementation requires PAE and pLDDT from structure predictor

        if pae_matrix is None or plddt_values is None:
            # Estimate from B-factors if PAE/pLDDT not available
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure("complex", pdb_path)

            b_factors = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            b_factors.append(atom.bfactor)

            avg_plddt = np.mean(b_factors) / 100.0
            pdockq = 0.7 * avg_plddt  # Rough approximation

        else:
            # Proper pDockQ calculation
            # This is a simplified version
            avg_plddt = np.mean(plddt_values) / 100.0
            avg_pae = np.mean(pae_matrix)
            pdockq = 0.724 / (1 + np.exp(0.052 * (avg_pae - 9.14)))

        return float(pdockq)