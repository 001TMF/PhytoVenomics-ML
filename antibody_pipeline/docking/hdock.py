"""
HDock integration for antibody-antigen docking

Based on: diffab-main/diffab/tools/dock/hdock.py
"""

import os
import shutil
import tempfile
import subprocess
import dataclasses as dc
from typing import List, Optional
from Bio import PDB
from Bio.PDB import Model as PDBModel

from ..preparation.renumber import renumber_antibody
from .base import DockingEngine


def fix_docked_pdb(pdb_path: str):
    """
    Fix HDock output PDB format issues.

    HDock sometimes produces malformed PDB files with incorrect
    line lengths. This function adds proper spacing.

    Args:
        pdb_path: Path to PDB file to fix
    """
    fixed = []
    with open(pdb_path, 'r') as f:
        for ln in f.readlines():
            if (ln.startswith('ATOM') or ln.startswith('HETATM')) and len(ln) == 56:
                fixed.append(ln[:-1] + ' 1.00  0.00              \n')
            else:
                fixed.append(ln)
    with open(pdb_path, 'w') as f:
        f.write(''.join(fixed))


class HDockEngine(DockingEngine):
    """
    HDock protein-protein docking engine.

    HDock is a fast docking algorithm that uses template-based
    and template-free modeling.

    Based on: diffab-main/diffab/tools/dock/hdock.py
    """

    def __init__(
        self,
        hdock_bin: str = './bin/hdock',
        createpl_bin: str = './bin/createpl',
    ):
        """
        Initialize HDock engine.

        Args:
            hdock_bin: Path to hdock executable
            createpl_bin: Path to createpl executable
        """
        super().__init__()
        self.hdock_bin = os.path.realpath(hdock_bin)
        self.createpl_bin = os.path.realpath(createpl_bin)
        self.tmpdir = tempfile.TemporaryDirectory()

        self._has_receptor = False
        self._has_ligand = False

        self._receptor_chains = []
        self._ligand_chains = []

    def __enter__(self):
        return self

    def __exit__(self, typ, value, traceback):
        self.tmpdir.cleanup()

    def set_receptor(self, pdb_path: str):
        """Set receptor (antigen) structure"""
        shutil.copyfile(pdb_path, os.path.join(self.tmpdir.name, 'receptor.pdb'))
        self._has_receptor = True

    def set_ligand(self, pdb_path: str):
        """Set ligand (antibody) structure"""
        shutil.copyfile(pdb_path, os.path.join(self.tmpdir.name, 'ligand.pdb'))
        self._has_ligand = True

    def _dump_complex_pdb(self) -> List[str]:
        """
        Combine receptor and docked ligand into complex PDB files.

        Returns:
            List of paths to complex PDB files
        """
        parser = PDB.PDBParser(QUIET=True)
        model_receptor = parser.get_structure(None, os.path.join(self.tmpdir.name, 'receptor.pdb'))[0]
        docked_pdb_path = os.path.join(self.tmpdir.name, 'ligand_docked.pdb')
        fix_docked_pdb(docked_pdb_path)
        structure_ligdocked = parser.get_structure(None, docked_pdb_path)

        pdb_io = PDB.PDBIO()
        paths = []
        for i, model_ligdocked in enumerate(structure_ligdocked):
            model_complex = PDBModel.Model(0)
            for chain in model_receptor:
                model_complex.add(chain.copy())
            for chain in model_ligdocked:
                model_complex.add(chain.copy())
            pdb_io.set_structure(model_complex)
            save_path = os.path.join(self.tmpdir.name, f"complex_{i}.pdb")
            pdb_io.save(save_path)
            paths.append(save_path)
        return paths

    def dock(self) -> List[str]:
        """
        Perform docking.

        Returns:
            List of paths to docked complex PDB files
        """
        if not (self._has_receptor and self._has_ligand):
            raise ValueError('Missing receptor or ligand.')
        subprocess.run(
            [self.hdock_bin, "receptor.pdb", "ligand.pdb"],
            cwd=self.tmpdir.name, check=True
        )
        subprocess.run(
            [self.createpl_bin, "Hdock.out", "ligand_docked.pdb"],
            cwd=self.tmpdir.name, check=True
        )
        return self._dump_complex_pdb()


@dc.dataclass
class DockSite:
    """
    Docking site specification.

    Attributes:
        chain: Chain identifier
        resseq: Residue sequence number
    """
    chain: str
    resseq: int


class HDockAntibody(HDockEngine):
    """
    HDock engine specialized for antibody-antigen docking.

    This variant:
    - Uses CDR H3 as the primary binding site
    - Supports epitope site constraints
    - Automatically renumbers antibody to Chothia scheme
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._heavy_chain_id = None
        self._epitope_sites: Optional[List[DockSite]] = None

    def set_ligand(self, pdb_path: str):
        """Disabled - use set_antibody instead"""
        raise NotImplementedError('Please use set_antibody')

    def set_receptor(self, pdb_path: str):
        """Disabled - use set_antigen instead"""
        raise NotImplementedError('Please use set_antigen')

    def set_antigen(self, pdb_path: str, epitope_sites: Optional[List[DockSite]] = None):
        """
        Set antigen structure with optional epitope constraints.

        Args:
            pdb_path: Path to antigen PDB file
            epitope_sites: Optional list of known epitope sites
        """
        super().set_receptor(pdb_path)
        self._epitope_sites = epitope_sites

    def set_antibody(self, pdb_path: str):
        """
        Set antibody structure.

        Automatically renumbers to Chothia scheme.

        Args:
            pdb_path: Path to antibody PDB file
        """
        heavy_chains, _ = renumber_antibody(
            pdb_path,
            os.path.join(self.tmpdir.name, 'ligand.pdb')
        )
        self._has_ligand = True
        self._heavy_chain_id = heavy_chains[0]

    def _prepare_lsite(self):
        """Prepare ligand (antibody) binding site file"""
        # Use CDR H3 (Chothia 95-102) as primary binding site
        lsite_content = f"95-102:{self._heavy_chain_id}\n"
        with open(os.path.join(self.tmpdir.name, 'lsite.txt'), 'w') as f:
            f.write(lsite_content)

    def _prepare_rsite(self):
        """Prepare receptor (antigen) binding site file"""
        rsite_content = ""
        for site in self._epitope_sites:
            rsite_content += f"{site.resseq}:{site.chain}\n"
        with open(os.path.join(self.tmpdir.name, 'rsite.txt'), 'w') as f:
            f.write(rsite_content)

    def dock(self) -> List[str]:
        """
        Perform antibody-antigen docking.

        Returns:
            List of paths to docked complex PDB files
        """
        if not (self._has_receptor and self._has_ligand):
            raise ValueError('Missing receptor or ligand.')

        self._prepare_lsite()

        # Build HDock command
        cmd_hdock = [self.hdock_bin, "receptor.pdb", "ligand.pdb", "-lsite", "lsite.txt"]
        if self._epitope_sites is not None:
            self._prepare_rsite()
            cmd_hdock += ["-rsite", "rsite.txt"]

        subprocess.run(
            cmd_hdock,
            cwd=self.tmpdir.name, check=True
        )

        # Build createpl command
        cmd_pl = [self.createpl_bin, "Hdock.out", "ligand_docked.pdb", "-lsite", "lsite.txt"]
        if self._epitope_sites is not None:
            cmd_pl += ["-rsite", "rsite.txt"]

        subprocess.run(
            cmd_pl,
            cwd=self.tmpdir.name, check=True
        )

        return self._dump_complex_pdb()