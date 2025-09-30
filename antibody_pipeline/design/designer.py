"""
RFdiffusion-based antibody design

Integrates RFdiffusion v3 for structure-based antibody design with
CDR masking and epitope conditioning.

Based on concepts from:
- diffab-main/diffab/tools/runner/design_for_pdb.py
- germinal-main/germinal/design/design.py
"""

import os
import enum
import subprocess
import tempfile
from typing import List, Optional, Dict, Tuple
from dataclasses import dataclass


class DesignMode(enum.Enum):
    """Antibody design modes"""
    SINGLE_CDR = "single_cdr"  # Design one CDR at a time
    MULTIPLE_CDRS = "multiple_cdrs"  # Design multiple CDRs simultaneously
    FULL = "full"  # Design entire antibody variable region
    CDR_H3_ONLY = "cdr_h3"  # Design only CDR H3 (most important for binding)


@dataclass
class DesignConfig:
    """Configuration for antibody design"""
    mode: DesignMode = DesignMode.CDR_H3_ONLY
    num_designs: int = 10
    num_steps: int = 50
    hotspot_res: Optional[List[str]] = None  # Epitope residues (e.g., ["A10", "A15"])
    temperature: float = 1.0
    noise_scale: float = 1.0


class AntibodyDesigner:
    """
    RFdiffusion-based antibody designer.

    This class provides an interface to RFdiffusion v3 for generating
    antibody structures conditioned on antigen epitopes.
    """

    def __init__(
        self,
        rfdiffusion_path: str,
        weights_path: Optional[str] = None,
        device: str = "cuda:0"
    ):
        """
        Initialize antibody designer.

        Args:
            rfdiffusion_path: Path to RFdiffusion installation
            weights_path: Path to model weights (optional)
            device: Device to use (cuda:0, cpu, etc.)
        """
        self.rfdiffusion_path = rfdiffusion_path
        self.weights_path = weights_path
        self.device = device

        # Validate RFdiffusion installation
        run_inference = os.path.join(rfdiffusion_path, "scripts", "run_inference.py")
        if not os.path.exists(run_inference):
            raise FileNotFoundError(
                f"RFdiffusion not found at {rfdiffusion_path}. "
                f"Expected {run_inference}"
            )
        self.run_inference = run_inference

    def design_antibody(
        self,
        antigen_pdb: str,
        output_dir: str,
        config: DesignConfig
    ) -> List[str]:
        """
        Design antibody structures for given antigen.

        Args:
            antigen_pdb: Path to antigen PDB file
            output_dir: Output directory for designs
            config: Design configuration

        Returns:
            List of paths to designed PDB files
        """
        os.makedirs(output_dir, exist_ok=True)

        # Build RFdiffusion command
        cmd = self._build_rfdiffusion_command(
            antigen_pdb,
            output_dir,
            config
        )

        # Run RFdiffusion
        print(f"Running RFdiffusion with command:\n{' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            cwd=self.rfdiffusion_path,
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"RFdiffusion failed with error:\n{result.stderr}"
            )

        # Collect output PDB files
        output_pdbs = []
        for fname in os.listdir(output_dir):
            if fname.endswith('.pdb'):
                output_pdbs.append(os.path.join(output_dir, fname))

        return sorted(output_pdbs)

    def _build_rfdiffusion_command(
        self,
        antigen_pdb: str,
        output_dir: str,
        config: DesignConfig
    ) -> List[str]:
        """
        Build RFdiffusion command line arguments.

        Args:
            antigen_pdb: Path to antigen PDB
            output_dir: Output directory
            config: Design configuration

        Returns:
            List of command arguments
        """
        cmd = [
            "python",
            self.run_inference.py,
            f"inference.output_prefix={output_dir}/design",
            f"inference.num_designs={config.num_designs}",
            f"diffuser.T={config.num_steps}",
            f"'contigmap.contigs=[A1-150/0 70-100]'",  # Antigen + antibody VH+VL
        ]

        # Add epitope hotspots if specified
        if config.hotspot_res:
            hotspot_str = ",".join(config.hotspot_res)
            cmd.append(f"'ppi.hotspot_res=[{hotspot_str}]'")

        # Add input PDB
        cmd.append(f"inference.input_pdb={antigen_pdb}")

        # Add weights path if specified
        if self.weights_path:
            cmd.append(f"inference.ckpt_override_path={self.weights_path}")

        return cmd

    def design_cdr_h3(
        self,
        antigen_pdb: str,
        framework_pdb: str,
        output_dir: str,
        num_designs: int = 10,
        epitope_residues: Optional[List[str]] = None
    ) -> List[str]:
        """
        Design CDR H3 region given antigen and antibody framework.

        CDR H3 is the most important region for antigen recognition
        and typically has the highest sequence diversity.

        Args:
            antigen_pdb: Path to antigen PDB
            framework_pdb: Path to antibody framework PDB
            output_dir: Output directory
            num_designs: Number of designs to generate
            epitope_residues: Known epitope residues (e.g., ["A10", "A15"])

        Returns:
            List of paths to designed PDB files
        """
        config = DesignConfig(
            mode=DesignMode.CDR_H3_ONLY,
            num_designs=num_designs,
            hotspot_res=epitope_residues
        )

        # For CDR H3 design, we need both antigen and framework
        # Create a temporary complex PDB
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            complex_pdb = tmp.name
            self._create_complex_pdb(antigen_pdb, framework_pdb, complex_pdb)

        try:
            output_pdbs = self.design_antibody(complex_pdb, output_dir, config)
        finally:
            os.unlink(complex_pdb)

        return output_pdbs

    def _create_complex_pdb(
        self,
        antigen_pdb: str,
        antibody_pdb: str,
        output_pdb: str
    ):
        """
        Create complex PDB from separate antigen and antibody files.

        Args:
            antigen_pdb: Path to antigen PDB
            antibody_pdb: Path to antibody PDB
            output_pdb: Path to output complex PDB
        """
        from Bio import PDB

        parser = PDB.PDBParser(QUIET=True)
        antigen_struct = parser.get_structure("antigen", antigen_pdb)
        antibody_struct = parser.get_structure("antibody", antibody_pdb)

        # Create new structure with both
        io = PDB.PDBIO()
        model = PDB.Model.Model(0)

        # Add antigen chains
        for chain in antigen_struct[0]:
            model.add(chain.copy())

        # Add antibody chains
        for chain in antibody_struct[0]:
            model.add(chain.copy())

        io.set_structure(model)
        io.save(output_pdb)

    def optimize_design(
        self,
        design_pdb: str,
        antigen_pdb: str,
        output_pdb: str,
        num_steps: int = 10
    ) -> str:
        """
        Optimize designed antibody using RFdiffusion partial diffusion.

        This performs additional refinement steps on an initial design.

        Args:
            design_pdb: Path to initial design PDB
            antigen_pdb: Path to antigen PDB
            output_pdb: Path to output optimized PDB
            num_steps: Number of optimization steps

        Returns:
            Path to optimized PDB
        """
        # Build optimization command
        cmd = [
            "python",
            self.run_inference,
            f"inference.input_pdb={design_pdb}",
            f"inference.output_prefix={os.path.splitext(output_pdb)[0]}",
            f"diffuser.T={num_steps}",
            f"diffuser.partial_T={num_steps}",  # Partial diffusion for refinement
        ]

        print(f"Optimizing design with {num_steps} steps...")
        result = subprocess.run(
            cmd,
            cwd=self.rfdiffusion_path,
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"RFdiffusion optimization failed:\n{result.stderr}"
            )

        return output_pdb


class DesignResult:
    """Container for design results with metadata"""

    def __init__(
        self,
        pdb_path: str,
        sequence: str,
        cdr_sequences: Dict[str, str],
        score: Optional[float] = None
    ):
        """
        Initialize design result.

        Args:
            pdb_path: Path to designed PDB
            sequence: Full antibody sequence
            cdr_sequences: Dictionary of CDR sequences (H1, H2, H3, L1, L2, L3)
            score: Design score (if available)
        """
        self.pdb_path = pdb_path
        self.sequence = sequence
        self.cdr_sequences = cdr_sequences
        self.score = score

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization"""
        return {
            'pdb_path': self.pdb_path,
            'sequence': self.sequence,
            'cdr_sequences': self.cdr_sequences,
            'score': self.score,
        }