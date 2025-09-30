"""
Antibody developability assessment

Based on: germinal-main/germinal/filters/filter_utils.py
"""

import numpy as np
from typing import Dict, List, Tuple
from Bio import PDB


class DevelopabilityFilter:
    """
    Assess antibody developability properties.

    Checks for:
    - Hydrophobic patches (SAP score)
    - Aggregation propensity
    - Post-translational modification sites
    - Charge patches
    - Framework mutations
    """

    def __init__(self):
        """Initialize developability filter"""
        pass

    def calculate_sap_score(
        self,
        pdb_path: str,
        chain_id: str,
        patch_radius: float = 9.0,
        avg_sasa_threshold: float = 15.0
    ) -> Dict[str, any]:
        """
        Calculate Spatial Aggregation Propensity (SAP) score.

        SAP identifies hydrophobic patches on the surface that
        may cause aggregation issues.

        Args:
            pdb_path: Path to PDB file
            chain_id: Chain to analyze
            patch_radius: Radius for patch detection (Angstroms)
            avg_sasa_threshold: SASA threshold for patch detection

        Returns:
            Dictionary with SAP score and patch information
        """
        from Bio.PDB import SASA

        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("antibody", pdb_path)

        # Calculate SASA
        sasa_calc = SASA.ShrakeRupley()
        sasa_calc.compute(structure, level="R")

        # Hydrophobicity scale (Kyte-Doolittle)
        hydrophobicity = {
            'ALA': 1.8, 'CYS': 2.5, 'ASP': -3.5, 'GLU': -3.5,
            'PHE': 2.8, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
            'LYS': -3.9, 'LEU': 3.8, 'MET': 1.9, 'ASN': -3.5,
            'PRO': -1.6, 'GLN': -3.5, 'ARG': -4.5, 'SER': -0.8,
            'THR': -0.7, 'VAL': 4.2, 'TRP': -0.9, 'TYR': -1.3
        }

        # Get residues from specified chain
        residues = []
        coords = []
        for model in structure:
            if chain_id in model:
                for residue in model[chain_id]:
                    if 'CA' in residue:
                        residues.append(residue)
                        coords.append(residue['CA'].coord)

        coords = np.array(coords)

        # Find hydrophobic patches
        patches = []
        for i, res in enumerate(residues):
            resname = res.get_resname()
            if resname not in hydrophobicity:
                continue

            # Check if residue is surface-exposed
            if res.sasa < avg_sasa_threshold:
                continue

            # Find neighbors within patch radius
            neighbors = []
            for j, other_res in enumerate(residues):
                if i != j:
                    dist = np.linalg.norm(coords[i] - coords[j])
                    if dist <= patch_radius:
                        neighbors.append(j)

            # Calculate average hydrophobicity of patch
            if neighbors:
                hydro_values = [hydrophobicity.get(res.get_resname(), 0.0)]
                for j in neighbors:
                    hydro_values.append(
                        hydrophobicity.get(residues[j].get_resname(), 0.0)
                    )
                avg_hydro = np.mean(hydro_values)

                if avg_hydro > 1.0:  # Hydrophobic patch threshold
                    patches.append({
                        'center': res.id[1],
                        'residue': resname,
                        'avg_hydrophobicity': avg_hydro,
                        'patch_size': len(neighbors) + 1
                    })

        # Calculate overall SAP score
        sap_score = sum(p['avg_hydrophobicity'] * p['patch_size'] for p in patches)

        return {
            'sap_score': sap_score,
            'num_patches': len(patches),
            'patches': patches,
        }

    def check_ptm_sites(
        self,
        sequence: str
    ) -> Dict[str, List[int]]:
        """
        Check for potential post-translational modification sites.

        Args:
            sequence: Amino acid sequence

        Returns:
            Dictionary of PTM sites
        """
        ptm_sites = {
            'N_glycosylation': [],  # N-X-S/T (X != P)
            'deamidation': [],  # N-G or N-S
            'isomerization': [],  # D-G or D-S
            'oxidation': [],  # M, W
        }

        # N-glycosylation: N-X-S/T where X is not P
        for i in range(len(sequence) - 2):
            if sequence[i] == 'N' and sequence[i + 1] != 'P' and sequence[i + 2] in ['S', 'T']:
                ptm_sites['N_glycosylation'].append(i)

        # Deamidation: N-G or N-S
        for i in range(len(sequence) - 1):
            if sequence[i] == 'N' and sequence[i + 1] in ['G', 'S']:
                ptm_sites['deamidation'].append(i)

        # Isomerization: D-G or D-S
        for i in range(len(sequence) - 1):
            if sequence[i] == 'D' and sequence[i + 1] in ['G', 'S']:
                ptm_sites['isomerization'].append(i)

        # Oxidation: M, W
        for i, aa in enumerate(sequence):
            if aa in ['M', 'W']:
                ptm_sites['oxidation'].append(i)

        return ptm_sites

    def check_charge_patches(
        self,
        pdb_path: str,
        chain_id: str,
        patch_radius: float = 10.0
    ) -> Dict[str, any]:
        """
        Identify charge patches that may affect binding or stability.

        Args:
            pdb_path: Path to PDB file
            chain_id: Chain to analyze
            patch_radius: Radius for patch detection (Angstroms)

        Returns:
            Dictionary with charge patch information
        """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("antibody", pdb_path)

        # Charge assignments
        positive = ['LYS', 'ARG', 'HIS']
        negative = ['ASP', 'GLU']

        # Get residues
        residues = []
        coords = []
        charges = []

        for model in structure:
            if chain_id in model:
                for residue in model[chain_id]:
                    if 'CA' in residue:
                        resname = residue.get_resname()
                        residues.append(residue)
                        coords.append(residue['CA'].coord)

                        if resname in positive:
                            charges.append(+1)
                        elif resname in negative:
                            charges.append(-1)
                        else:
                            charges.append(0)

        coords = np.array(coords)
        charges = np.array(charges)

        # Find charge patches
        patches = []
        for i, res in enumerate(residues):
            if charges[i] == 0:
                continue

            # Find neighbors within patch radius
            neighbors = []
            for j in range(len(residues)):
                if i != j:
                    dist = np.linalg.norm(coords[i] - coords[j])
                    if dist <= patch_radius:
                        neighbors.append(j)

            # Calculate net charge of patch
            if neighbors:
                patch_charges = [charges[i]] + [charges[j] for j in neighbors]
                net_charge = sum(patch_charges)
                abs_charge = abs(net_charge)

                if abs_charge >= 3:  # Significant charge patch
                    patches.append({
                        'center': res.id[1],
                        'residue': res.get_resname(),
                        'net_charge': net_charge,
                        'patch_size': len(neighbors) + 1
                    })

        return {
            'num_charge_patches': len(patches),
            'charge_patches': patches,
        }

    def count_framework_mutations(
        self,
        sequence: str,
        reference_sequence: str,
        cdr_positions: List[int]
    ) -> Dict[str, any]:
        """
        Count mutations in framework regions.

        Args:
            sequence: Designed sequence
            reference_sequence: Reference (template) sequence
            cdr_positions: Positions of CDR residues

        Returns:
            Dictionary with mutation information
        """
        if len(sequence) != len(reference_sequence):
            raise ValueError("Sequences must have same length")

        mutations = []
        framework_mutations = []

        for i, (aa1, aa2) in enumerate(zip(sequence, reference_sequence)):
            if aa1 != aa2:
                mutation = {
                    'position': i,
                    'from': aa2,
                    'to': aa1,
                    'is_framework': i not in cdr_positions
                }
                mutations.append(mutation)

                if mutation['is_framework']:
                    framework_mutations.append(mutation)

        return {
            'num_mutations': len(mutations),
            'num_framework_mutations': len(framework_mutations),
            'mutations': mutations,
            'framework_mutations': framework_mutations,
        }

    def assess_developability(
        self,
        pdb_path: str,
        sequence: str,
        chain_id: str,
        reference_sequence: Optional[str] = None,
        cdr_positions: Optional[List[int]] = None
    ) -> Dict[str, any]:
        """
        Comprehensive developability assessment.

        Args:
            pdb_path: Path to PDB file
            sequence: Amino acid sequence
            chain_id: Chain to analyze
            reference_sequence: Reference sequence for mutation analysis
            cdr_positions: CDR residue positions

        Returns:
            Dictionary with all developability metrics
        """
        results = {}

        # SAP score
        results['sap'] = self.calculate_sap_score(pdb_path, chain_id)

        # PTM sites
        results['ptm_sites'] = self.check_ptm_sites(sequence)

        # Charge patches
        results['charge_patches'] = self.check_charge_patches(pdb_path, chain_id)

        # Framework mutations
        if reference_sequence and cdr_positions:
            results['mutations'] = self.count_framework_mutations(
                sequence, reference_sequence, cdr_positions
            )

        # Overall assessment
        issues = []
        if results['sap']['sap_score'] > 50:
            issues.append('High SAP score (aggregation risk)')
        if len(results['ptm_sites']['N_glycosylation']) > 0:
            issues.append('N-glycosylation sites present')
        if len(results['ptm_sites']['deamidation']) > 2:
            issues.append('Multiple deamidation sites')
        if results['charge_patches']['num_charge_patches'] > 3:
            issues.append('Multiple charge patches')

        results['developable'] = len(issues) == 0
        results['issues'] = issues

        return results