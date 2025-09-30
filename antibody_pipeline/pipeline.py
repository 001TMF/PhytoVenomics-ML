"""
Main antibody design pipeline orchestrator

Integrates all modules into end-to-end workflow:
1. Preparation (renumbering, CDR identification)
2. Docking (antibody-antigen complex generation)
3. Design (RFdiffusion structure generation)
4. Optimization (IgLM sequence refinement)
5. Filtering (quality assessment and ranking)
"""

import os
import json
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict

from .preparation.renumber import AntibodyRenumberer
from .preparation.parser import preprocess_antibody_structure
from .preparation.cdr_identifier import get_cdr_info, extract_cdr_sequences
from .docking.hdock import HDockAntibody, DockSite
from .design.designer import AntibodyDesigner, DesignConfig, DesignMode
from .optimization.iglm_optimizer import IgLMOptimizer
from .filtering.structure_filters import StructureFilter
from .filtering.interface_metrics import InterfaceMetrics
from .filtering.developability import DevelopabilityFilter


@dataclass
class PipelineConfig:
    """Configuration for the antibody design pipeline"""

    # Input files - THREE MODES:
    # Mode 1: Separate antigen + framework (will dock)
    antigen_pdb: Optional[str] = None
    framework_pdb: Optional[str] = None

    # Mode 2: Pre-docked complex (skip docking)
    complex_pdb: Optional[str] = None
    antibody_chains: Optional[List[str]] = None  # e.g., ["H", "L"]
    antigen_chains: Optional[List[str]] = None   # e.g., ["A"]

    # Mode 3: Antigen only (de novo design)
    # Just provide antigen_pdb and no framework_pdb or complex_pdb

    # Output directory
    output_dir: str = "./results"

    # Design settings
    design_mode: str = "cdr_h3"  # single_cdr, multiple_cdrs, full, cdr_h3
    num_designs: int = 10
    num_design_steps: int = 50

    # Epitope information
    epitope_residues: Optional[List[str]] = None  # e.g., ["A10", "A15"]

    # Docking settings
    use_docking: bool = True
    hdock_bin: str = "./bin/hdock"
    createpl_bin: str = "./bin/createpl"

    # RFdiffusion settings
    rfdiffusion_path: str = "./RFdiffusion"
    rfdiffusion_weights: Optional[str] = None

    # IgLM optimization
    use_iglm_optimization: bool = True
    iglm_species: str = "[HUMAN]"
    iglm_temperature: float = 1.0

    # Filtering thresholds
    max_clashes: int = 100
    min_plddt: float = 70.0
    min_pdockq: float = 0.23
    max_sap_score: float = 50.0
    min_cdr_interface_pct: float = 60.0

    # Structure prediction
    structure_predictor: str = "chai"  # "chai" or "af3"
    chai_path: Optional[str] = None
    af3_path: Optional[str] = None

    # Device
    device: str = "cuda:0"


@dataclass
class DesignCandidate:
    """Container for design candidate with all metrics"""
    design_id: str
    sequence: str
    pdb_path: str
    cdr_sequences: Dict[str, str]

    # IgLM scores
    iglm_log_likelihood: Optional[float] = None
    iglm_perplexity: Optional[float] = None

    # Structure metrics
    plddt: Optional[float] = None
    ptm: Optional[float] = None
    iptm: Optional[float] = None
    num_clashes: Optional[int] = None

    # Interface metrics
    num_contacts: Optional[int] = None
    bsa: Optional[float] = None
    pdockq: Optional[float] = None
    cdr_interface_pct: Optional[float] = None

    # Developability
    sap_score: Optional[float] = None
    num_ptm_sites: Optional[int] = None
    num_charge_patches: Optional[int] = None

    # Overall
    passed_filters: bool = False
    rank: Optional[int] = None


class AntibodyDesignPipeline:
    """
    Complete end-to-end antibody design pipeline.

    Combines DiffAb preparation, HDock docking, RFdiffusion design,
    IgLM optimization, and Germinal-style filtering.
    """

    def __init__(self, config: PipelineConfig):
        """
        Initialize pipeline.

        Args:
            config: Pipeline configuration
        """
        self.config = config
        os.makedirs(config.output_dir, exist_ok=True)

        # Initialize modules
        self.renumberer = AntibodyRenumberer()

        if config.use_docking:
            self.docker = None  # Initialized on demand

        self.designer = AntibodyDesigner(
            rfdiffusion_path=config.rfdiffusion_path,
            weights_path=config.rfdiffusion_weights,
            device=config.device
        )

        if config.use_iglm_optimization:
            self.iglm = IgLMOptimizer(
                species=config.iglm_species,
                temperature=config.iglm_temperature,
                device=config.device
            )

        self.structure_filter = StructureFilter(
            predictor=config.structure_predictor,
            chai_path=config.chai_path,
            af3_path=config.af3_path
        )

        self.interface_metrics = InterfaceMetrics()
        self.developability_filter = DevelopabilityFilter()

        self.candidates: List[DesignCandidate] = []

    def run(self) -> List[DesignCandidate]:
        """
        Run complete pipeline.

        Returns:
            List of ranked design candidates
        """
        print("=" * 80)
        print("ANTIBODY DESIGN PIPELINE")
        print("=" * 80)

        # Step 1: Preparation
        print("\n[1/6] Preparing inputs...")
        prepared_data = self._prepare_inputs()

        # Step 2: Docking (if needed)
        if prepared_data['mode'] == 'pre-docked':
            print("\n[2/6] Using pre-docked complex (skipping docking)...")
            docked_pdbs = [prepared_data['complex_pdb']]
        elif prepared_data['mode'] == 'separate' and self.config.use_docking:
            print("\n[2/6] Docking antibody-antigen complex...")
            docked_pdbs = self._dock_complex(prepared_data)
        else:
            print("\n[2/6] Skipping docking...")
            docked_pdbs = [prepared_data['antigen_pdb']]

        # Step 3: Design
        print(f"\n[3/6] Generating {self.config.num_designs} designs...")
        design_pdbs = self._generate_designs(docked_pdbs[0])

        # Step 4: Optimization
        if self.config.use_iglm_optimization:
            print("\n[4/6] Optimizing sequences with IgLM...")
            optimized_designs = self._optimize_sequences(design_pdbs)
        else:
            print("\n[4/6] Skipping IgLM optimization...")
            optimized_designs = design_pdbs

        # Step 5: Filtering
        print("\n[5/6] Filtering and scoring designs...")
        self._filter_designs(optimized_designs)

        # Step 6: Ranking
        print("\n[6/6] Ranking designs...")
        ranked_candidates = self._rank_designs()

        # Save results
        self._save_results(ranked_candidates)

        print("\n" + "=" * 80)
        print(f"Pipeline complete! {len(ranked_candidates)} designs passed filters.")
        print(f"Results saved to: {self.config.output_dir}")
        print("=" * 80)

        return ranked_candidates

    def _prepare_inputs(self) -> Dict:
        """Prepare and validate input structures"""
        antigen_dir = os.path.join(self.config.output_dir, "prepared")
        os.makedirs(antigen_dir, exist_ok=True)

        # Determine input mode
        if self.config.complex_pdb:
            # MODE 2: Pre-docked complex
            print("  Using pre-docked complex")

            # Split complex into antibody and antigen
            antibody_pdb, antigen_pdb = self._split_complex(
                self.config.complex_pdb,
                self.config.antibody_chains,
                self.config.antigen_chains,
                antigen_dir
            )

            # Renumber antibody
            framework_renumbered = os.path.join(antigen_dir, "antibody_chothia.pdb")
            result = self.renumberer.renumber(antibody_pdb, framework_renumbered)
            print(f"  Antibody renumbered: {result['heavy_chains']} (H), {result['light_chains']} (L)")

            return {
                'antigen_pdb': antigen_pdb,
                'framework_pdb': framework_renumbered,
                'complex_pdb': self.config.complex_pdb,
                'mode': 'pre-docked'
            }

        elif self.config.framework_pdb:
            # MODE 1: Separate antigen + framework (will dock)
            print("  Using separate antigen and framework")
            framework_renumbered = os.path.join(antigen_dir, "framework_chothia.pdb")
            result = self.renumberer.renumber(
                self.config.framework_pdb,
                framework_renumbered
            )
            print(f"  Framework renumbered: {result['heavy_chains']} (H), {result['light_chains']} (L)")

            return {
                'antigen_pdb': self.config.antigen_pdb,
                'framework_pdb': framework_renumbered,
                'mode': 'separate'
            }

        else:
            # MODE 3: Antigen only (de novo)
            print("  De novo design mode (antigen only)")
            return {
                'antigen_pdb': self.config.antigen_pdb,
                'framework_pdb': None,
                'mode': 'de-novo'
            }

    def _split_complex(
        self,
        complex_pdb: str,
        antibody_chains: List[str],
        antigen_chains: List[str],
        output_dir: str
    ) -> Tuple[str, str]:
        """
        Split pre-docked complex into separate antibody and antigen PDB files.

        Args:
            complex_pdb: Path to complex PDB
            antibody_chains: List of antibody chain IDs (e.g., ["H", "L"])
            antigen_chains: List of antigen chain IDs (e.g., ["A"])
            output_dir: Output directory

        Returns:
            Tuple of (antibody_pdb_path, antigen_pdb_path)
        """
        from Bio import PDB

        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("complex", complex_pdb)

        # Create separate structures
        antibody_structure = PDB.Structure.Structure("antibody")
        antigen_structure = PDB.Structure.Structure("antigen")

        for model in structure:
            ab_model = PDB.Model.Model(0)
            ag_model = PDB.Model.Model(0)

            for chain in model:
                if chain.id in antibody_chains:
                    ab_model.add(chain.copy())
                elif chain.id in antigen_chains:
                    ag_model.add(chain.copy())

            if len(ab_model):
                antibody_structure.add(ab_model)
            if len(ag_model):
                antigen_structure.add(ag_model)

        # Save separate PDB files
        io = PDB.PDBIO()

        antibody_pdb = os.path.join(output_dir, "antibody_from_complex.pdb")
        io.set_structure(antibody_structure)
        io.save(antibody_pdb)

        antigen_pdb = os.path.join(output_dir, "antigen_from_complex.pdb")
        io.set_structure(antigen_structure)
        io.save(antigen_pdb)

        print(f"  Split complex: antibody chains {antibody_chains}, antigen chains {antigen_chains}")

        return antibody_pdb, antigen_pdb

    def _dock_complex(self, prepared_data: Dict) -> List[str]:
        """Dock antibody-antigen complex using HDock"""
        # Convert epitope residues to DockSite format
        epitope_sites = None
        if self.config.epitope_residues:
            epitope_sites = []
            for res_str in self.config.epitope_residues:
                chain = res_str[0]
                resseq = int(res_str[1:])
                epitope_sites.append(DockSite(chain=chain, resseq=resseq))

        with HDockAntibody(
            hdock_bin=self.config.hdock_bin,
            createpl_bin=self.config.createpl_bin
        ) as docker:
            docker.set_antigen(prepared_data['antigen_pdb'], epitope_sites=epitope_sites)
            docker.set_antibody(prepared_data['framework_pdb'])
            docked_complexes = docker.dock()

        # Save docked complexes
        dock_dir = os.path.join(self.config.output_dir, "docked")
        os.makedirs(dock_dir, exist_ok=True)

        saved_paths = []
        for i, complex_pdb in enumerate(docked_complexes[:3]):  # Keep top 3
            dest_path = os.path.join(dock_dir, f"docked_{i}.pdb")
            import shutil
            shutil.copy(complex_pdb, dest_path)
            saved_paths.append(dest_path)
            print(f"  Saved docked complex: {dest_path}")

        return saved_paths

    def _generate_designs(self, template_pdb: str) -> List[str]:
        """Generate antibody designs using RFdiffusion"""
        design_dir = os.path.join(self.config.output_dir, "designs")
        os.makedirs(design_dir, exist_ok=True)

        # Configure design
        design_config = DesignConfig(
            mode=DesignMode(self.config.design_mode),
            num_designs=self.config.num_designs,
            num_steps=self.config.num_design_steps,
            hotspot_res=self.config.epitope_residues
        )

        # Generate designs
        if self.config.framework_pdb:
            # Design with framework
            design_pdbs = self.designer.design_cdr_h3(
                antigen_pdb=self.config.antigen_pdb,
                framework_pdb=self.config.framework_pdb,
                output_dir=design_dir,
                num_designs=self.config.num_designs,
                epitope_residues=self.config.epitope_residues
            )
        else:
            # De novo design
            design_pdbs = self.designer.design_antibody(
                antigen_pdb=self.config.antigen_pdb,
                output_dir=design_dir,
                config=design_config
            )

        print(f"  Generated {len(design_pdbs)} designs")
        return design_pdbs

    def _optimize_sequences(self, design_pdbs: List[str]) -> List[str]:
        """Optimize sequences using IgLM"""
        from Bio import PDB, SeqIO

        optimized_pdbs = []
        for i, pdb_path in enumerate(design_pdbs):
            # Extract sequence
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure("design", pdb_path)

            # Get heavy chain sequence (assuming chain H)
            sequence = ""
            for model in structure:
                if 'H' in model:
                    for residue in model['H']:
                        if PDB.is_aa(residue):
                            sequence += PDB.Polypeptide.three_to_one(residue.get_resname())

            # Score with IgLM
            score = self.iglm.score_sequence(sequence, chain="[HEAVY]")
            print(f"  Design {i}: IgLM LL = {score['log_likelihood']:.2f}, PPL = {score['perplexity']:.2f}")

            optimized_pdbs.append(pdb_path)

        return optimized_pdbs

    def _filter_designs(self, design_pdbs: List[str]):
        """Filter and score all designs"""
        from Bio import PDB

        for i, pdb_path in enumerate(design_pdbs):
            design_id = f"design_{i:03d}"
            print(f"\n  Evaluating {design_id}...")

            try:
                # Extract sequence
                parser = PDB.PDBParser(QUIET=True)
                structure = parser.get_structure("design", pdb_path)

                sequence = ""
                for model in structure:
                    if 'H' in model:
                        for residue in model['H']:
                            if PDB.is_aa(residue):
                                sequence += PDB.Polypeptide.three_to_one(residue.get_resname())

                # Create candidate
                candidate = DesignCandidate(
                    design_id=design_id,
                    sequence=sequence,
                    pdb_path=pdb_path,
                    cdr_sequences={}  # Would extract properly in full implementation
                )

                # IgLM score
                if self.config.use_iglm_optimization:
                    iglm_score = self.iglm.score_sequence(sequence, chain="[HEAVY]")
                    candidate.iglm_log_likelihood = iglm_score['log_likelihood']
                    candidate.iglm_perplexity = iglm_score['perplexity']

                # Structure validation
                passed, metrics = self.structure_filter.validate_structure(
                    pdb_path,
                    max_clashes=self.config.max_clashes,
                    min_plddt=self.config.min_plddt
                )
                candidate.num_clashes = metrics['num_clashes']
                candidate.plddt = metrics['avg_plddt']

                # Developability
                dev_results = self.developability_filter.assess_developability(
                    pdb_path, sequence, chain_id='H'
                )
                candidate.sap_score = dev_results['sap']['sap_score']
                candidate.passed_filters = dev_results['developable'] and passed

                self.candidates.append(candidate)

                print(f"    Clashes: {candidate.num_clashes}, pLDDT: {candidate.plddt:.1f}")
                print(f"    SAP: {candidate.sap_score:.1f}, Passed: {candidate.passed_filters}")

            except Exception as e:
                print(f"    Error processing {design_id}: {e}")
                continue

    def _rank_designs(self) -> List[DesignCandidate]:
        """Rank designs by composite score"""
        # Filter passing candidates
        passing = [c for c in self.candidates if c.passed_filters]

        # Sort by composite score
        # Higher IgLM LL, higher pLDDT, lower SAP
        for candidate in passing:
            score = 0.0
            if candidate.iglm_log_likelihood:
                score += candidate.iglm_log_likelihood
            if candidate.plddt:
                score += candidate.plddt / 100.0
            if candidate.sap_score:
                score -= candidate.sap_score / 50.0
            candidate.rank_score = score

        passing.sort(key=lambda c: c.rank_score, reverse=True)

        # Assign ranks
        for i, candidate in enumerate(passing):
            candidate.rank = i + 1

        return passing

    def _save_results(self, candidates: List[DesignCandidate]):
        """Save results to JSON"""
        results = {
            'config': asdict(self.config),
            'num_designs': len(self.candidates),
            'num_passing': len(candidates),
            'candidates': [asdict(c) for c in candidates]
        }

        results_path = os.path.join(self.config.output_dir, "results.json")
        with open(results_path, 'w') as f:
            json.dump(results, f, indent=2)

        print(f"\n  Results saved: {results_path}")