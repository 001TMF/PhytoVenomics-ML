#!/usr/bin/env python
# antibody_design/antibody_generator.py

import logging
import numpy as np
import torch
from typing import Dict, List, Union, Tuple, Optional
import os
import random
import json
from pathlib import Path

logger = logging.getLogger("Phytovenomics.HumanAntibodyDesigner")

class HumanAntibodyDesigner:
    """
    Designs human antibodies against identified toxin epitopes using protein language models.
    """
    def __init__(self, config: Dict):
        """
        Initialize the HumanAntibodyDesigner with configuration settings.
        
        Args:
            config: Dictionary containing configuration parameters
        """
        self.config = config
        
        # Load model
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model_name = self.config.get("model_name", "esmif_v1")
        
        # Define CDR regions
        self.cdr_definitions = self.config.get("cdr_definitions", [
            {"name": "CDR-H1", "start": 26, "end": 32},
            {"name": "CDR-H2", "start": 52, "end": 56},
            {"name": "CDR-H3", "start": 95, "end": 102},
            {"name": "CDR-L1", "start": 24, "end": 34},
            {"name": "CDR-L2", "start": 50, "end": 56},
            {"name": "CDR-L3", "start": 89, "end": 97},
        ])
        
        # Load template framework
        framework_path = self.config.get("framework_template", "human_germline")
        self.framework = self._load_framework(framework_path)
        
        # Try to load models if available
        try:
            if torch.cuda.is_available():
                logger.info("Using GPU for antibody design")
            else:
                logger.info("Using CPU for antibody design")
                
            logger.info(f"Loading antibody design model {self.model_name}...")
            # In a production implementation, actual model loading would happen here
            logger.info("Model loaded successfully (placeholder)")
        except Exception as e:
            logger.error(f"Error loading model: {e}")
            
    def _load_framework(self, framework_path: str) -> Dict:
        """Load antibody framework template"""
        try:
            # Check if it's one of the built-in templates
            built_in = {
                "human_germline": {
                    "heavy": {
                        "sequence": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKGGGGSYAMDYWGQGTLVTVSS",
                        "name": "IGHV3-23*01",
                        "CDRs": {
                            "CDR-H1": (26, 32),
                            "CDR-H2": (52, 56),
                            "CDR-H3": (95, 102)
                        }
                    },
                    "light": {
                        "sequence": "DIQMTQSPSSLSASVGDRVTITCRASQSVSSAVAWYQQKPGKAPKLLIYSASSLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQYSGSPLTFGGGTKVEIK",
                        "name": "IGKV1-39*01",
                        "CDRs": {
                            "CDR-L1": (24, 34),
                            "CDR-L2": (50, 56),
                            "CDR-L3": (89, 97)
                        }
                    }
                }
            }
            
            if framework_path in built_in:
                logger.info(f"Using built-in framework template: {framework_path}")
                return built_in[framework_path]
            
            # Otherwise, try to load from a file
            if os.path.exists(framework_path):
                with open(framework_path, 'r') as f:
                    data = json.load(f)
                return data
            
            logger.warning(f"Framework path {framework_path} not found, using default human framework")
            return built_in["human_germline"]
            
        except Exception as e:
            logger.error(f"Error loading framework template: {e}")
            # Return a minimal default framework
            return {
                "heavy": {
                    "sequence": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKGGGGSYAMDYWGQGTLVTVSS",
                    "name": "IGHV3-23*01"
                },
                "light": {
                    "sequence": "DIQMTQSPSSLSASVGDRVTITCRASQSVSSAVAWYQQKPGKAPKLLIYSASSLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQYSGSPLTFGGGTKVEIK",
                    "name": "IGKV1-39*01"
                }
            }
    
    def generate_bnab(self, target: Dict) -> Dict:
        """
        Design a broadly neutralizing antibody against a conserved epitope.
        
        Args:
            target: Dictionary containing epitope information
            
        Returns:
            Dictionary containing the designed antibody information
        """
        logger.info(f"Generating broadly neutralizing antibody for target {target['id']}")
        
        try:
            # Create a unique ID for this antibody
            antibody_id = f"bnab_{target['id']}_{random.randint(1000, 9999)}"
            
            # Get the epitope sequence
            epitope_seq = target["sequence"]
            
            # Generate CDR sequences based on the epitope
            cdr_sequences = self._design_cdrs_for_epitope(epitope_seq, is_bnab=True)
            
            # Create the full antibody by grafting CDRs onto the framework
            heavy_chain, light_chain = self._graft_cdrs_to_framework(cdr_sequences)
            
            # Calculate binding energy (prediction)
            binding_energy = self._predict_binding_energy(heavy_chain, light_chain, epitope_seq)
            
            # Calculate developability score
            developability = self._calculate_developability(heavy_chain, light_chain)
            
            # Create antibody object
            antibody = {
                "id": antibody_id,
                "type": "broadly_neutralizing",
                "target": target["id"],
                "heavy_chain": {
                    "sequence": heavy_chain,
                    "framework": self.framework["heavy"]["name"],
                    "CDRs": {
                        cdr["name"]: cdr_sequences[cdr["name"]]
                        for cdr in self.cdr_definitions
                        if cdr["name"].startswith("CDR-H")
                    }
                },
                "light_chain": {
                    "sequence": light_chain,
                    "framework": self.framework["light"]["name"],
                    "CDRs": {
                        cdr["name"]: cdr_sequences[cdr["name"]]
                        for cdr in self.cdr_definitions
                        if cdr["name"].startswith("CDR-L")
                    }
                },
                "predicted_binding": {
                    "energy": binding_energy,
                    "affinity_nm": self._convert_energy_to_affinity(binding_energy)
                },
                "developability": developability
            }
            
            return antibody
            
        except Exception as e:
            logger.exception(f"Error generating broadly neutralizing antibody: {e}")
            # Return a minimal placeholder antibody
            return {
                "id": f"error_bnab_{target['id']}",
                "type": "broadly_neutralizing",
                "target": target["id"],
                "heavy_chain": {"sequence": self.framework["heavy"]["sequence"]},
                "light_chain": {"sequence": self.framework["light"]["sequence"]},
                "error": str(e)
            }
    
    def generate_specific_ab(self, target: Dict) -> Dict:
        """
        Design a specific antibody against a unique epitope on a toxin.
        
        Args:
            target: Dictionary containing epitope information
            
        Returns:
            Dictionary containing the designed antibody information
        """
        logger.info(f"Generating specific antibody for target {target['id']}")
        
        try:
            # Create a unique ID for this antibody
            antibody_id = f"specific_{target['id']}_{random.randint(1000, 9999)}"
            
            # Get the epitope sequence
            epitope_seq = target["sequence"]
            
            # Generate CDR sequences optimized for specificity
            cdr_sequences = self._design_cdrs_for_epitope(epitope_seq, is_bnab=False)
            
            # Create the full antibody by grafting CDRs onto the framework
            heavy_chain, light_chain = self._graft_cdrs_to_framework(cdr_sequences)
            
            # Calculate binding energy (prediction)
            binding_energy = self._predict_binding_energy(heavy_chain, light_chain, epitope_seq)
            
            # Calculate developability score
            developability = self._calculate_developability(heavy_chain, light_chain)
            
            # Calculate specificity score
            specificity = target.get("uniqueness_score", 0.8)
            
            # Create antibody object
            antibody = {
                "id": antibody_id,
                "type": "specific",
                "target": target["id"],
                "heavy_chain": {
                    "sequence": heavy_chain,
                    "framework": self.framework["heavy"]["name"],
                    "CDRs": {
                        cdr["name"]: cdr_sequences[cdr["name"]]
                        for cdr in self.cdr_definitions
                        if cdr["name"].startswith("CDR-H")
                    }
                },
                "light_chain": {
                    "sequence": light_chain,
                    "framework": self.framework["light"]["name"],
                    "CDRs": {
                        cdr["name"]: cdr_sequences[cdr["name"]]
                        for cdr in self.cdr_definitions
                        if cdr["name"].startswith("CDR-L")
                    }
                },
                "predicted_binding": {
                    "energy": binding_energy,
                    "affinity_nm": self._convert_energy_to_affinity(binding_energy)
                },
                "developability": developability,
                "specificity": specificity
            }
            
            return antibody
            
        except Exception as e:
            logger.exception(f"Error generating specific antibody: {e}")
            # Return a minimal placeholder antibody
            return {
                "id": f"error_specific_{target['id']}",
                "type": "specific",
                "target": target["id"],
                "heavy_chain": {"sequence": self.framework["heavy"]["sequence"]},
                "light_chain": {"sequence": self.framework["light"]["sequence"]},
                "error": str(e)
            }
            
    def cdr_optimization(self, antibody: Dict, target: Dict) -> Dict:
        """
        Optimize CDR sequences to improve binding while maintaining developability.
        
        Args:
            antibody: Dictionary containing antibody information
            target: Dictionary containing epitope information
            
        Returns:
            Updated antibody dictionary with optimized sequences
        """
        logger.info(f"Optimizing CDRs for antibody {antibody['id']}")
        
        try:
            # Extract epitope sequence
            epitope_seq = target["sequence"]
            
            # Extract current CDR sequences
            current_cdrs = {}
            for cdr in self.cdr_definitions:
                if cdr["name"].startswith("CDR-H") and "heavy_chain" in antibody and "CDRs" in antibody["heavy_chain"]:
                    current_cdrs[cdr["name"]] = antibody["heavy_chain"]["CDRs"].get(cdr["name"])
                elif cdr["name"].startswith("CDR-L") and "light_chain" in antibody and "CDRs" in antibody["light_chain"]:
                    current_cdrs[cdr["name"]] = antibody["light_chain"]["CDRs"].get(cdr["name"])
            
            # If we don't have CDRs, extract them from the sequences
            if not current_cdrs or None in current_cdrs.values():
                current_cdrs = self._extract_cdrs_from_sequences(
                    antibody["heavy_chain"]["sequence"],
                    antibody["light_chain"]["sequence"]
                )
            
            # Make focused optimizations to CDR-H3 (most important for binding)
            if "CDR-H3" in current_cdrs and current_cdrs["CDR-H3"]:
                # Generate candidates for CDR-H3
                candidates = self._generate_cdr_variants(current_cdrs["CDR-H3"], epitope_seq, "CDR-H3")
                
                # Score candidates
                best_candidate = current_cdrs["CDR-H3"]
                best_score = self._score_cdr_epitope_binding(current_cdrs["CDR-H3"], epitope_seq)
                
                for candidate in candidates:
                    score = self._score_cdr_epitope_binding(candidate, epitope_seq)
                    if score > best_score:
                        best_score = score
                        best_candidate = candidate
                
                # Update CDR-H3
                current_cdrs["CDR-H3"] = best_candidate
            
            # Regenerate full antibody sequences
            heavy_chain, light_chain = self._graft_cdrs_to_framework(current_cdrs)
            
            # Update binding energy
            binding_energy = self._predict_binding_energy(heavy_chain, light_chain, epitope_seq)
            
            # Update developability
            developability = self._calculate_developability(heavy_chain, light_chain)
            
            # Update antibody object
            optimized_antibody = antibody.copy()
            optimized_antibody["heavy_chain"] = {
                "sequence": heavy_chain,
                "framework": antibody["heavy_chain"]["framework"] if "framework" in antibody["heavy_chain"] else self.framework["heavy"]["name"],
                "CDRs": {
                    cdr["name"]: current_cdrs[cdr["name"]]
                    for cdr in self.cdr_definitions
                    if cdr["name"].startswith("CDR-H") and cdr["name"] in current_cdrs
                }
            }
            optimized_antibody["light_chain"] = {
                "sequence": light_chain,
                "framework": antibody["light_chain"]["framework"] if "framework" in antibody["light_chain"] else self.framework["light"]["name"],
                "CDRs": {
                    cdr["name"]: current_cdrs[cdr["name"]]
                    for cdr in self.cdr_definitions
                    if cdr["name"].startswith("CDR-L") and cdr["name"] in current_cdrs
                }
            }
            optimized_antibody["predicted_binding"] = {
                "energy": binding_energy,
                "affinity_nm": self._convert_energy_to_affinity(binding_energy)
            }
            optimized_antibody["developability"] = developability
            optimized_antibody["optimized"] = True
            
            return optimized_antibody
            
        except Exception as e:
            logger.exception(f"Error optimizing CDRs: {e}")
            # Return the original antibody
            antibody["error_optimization"] = str(e)
            return antibody

    def _design_cdrs_for_epitope(self, epitope_seq: str, is_bnab: bool = False) -> Dict[str, str]:
        """Generate CDR sequences for the given epitope"""
        # This is a simplified implementation that would be replaced with ML-based design
        # For now, we'll use template-based design with customizations
        
        # Define amino acid properties
        charged_pos = "KR"  # Positively charged
        charged_neg = "DE"  # Negatively charged
        polar = "NQST"      # Polar
        hydrophobic = "AILMFWYV"  # Hydrophobic
        h_bond = "YST"      # H-bond formers
        aromatic = "FYW"    # Aromatic
        
        # Analyze epitope composition
        epi_len = len(epitope_seq)
        charged_pos_count = sum(1 for aa in epitope_seq if aa in charged_pos)
        charged_neg_count = sum(1 for aa in epitope_seq if aa in charged_neg)
        polar_count = sum(1 for aa in epitope_seq if aa in polar)
        hydrophobic_count = sum(1 for aa in epitope_seq if aa in hydrophobic)
        
        # Define template CDRs
        template_cdrs = {
            "CDR-H1": "GFTFSSYAMS",
            "CDR-H2": "AISGSGGSTYYADSVKG",
            "CDR-H3": "AKGGRRTYYYGSYAMDY" if is_bnab else "AKGGGGSYAMDY",
            "CDR-L1": "RASQSVSSAVA",
            "CDR-L2": "SASSLYS", 
            "CDR-L3": "QQYSGSPL"
        }
        
        # Customize CDR-H3 (most important for binding)
        cdr_h3 = list(template_cdrs["CDR-H3"])
        
        # For negatively charged epitopes, add positive charges
        if charged_neg_count / epi_len > 0.2:
            for i in range(3, min(len(cdr_h3) - 4, 10)):
                if random.random() < 0.4:  # 40% chance to modify
                    cdr_h3[i] = random.choice(charged_pos)
        
        # For positively charged epitopes, add negative charges
        elif charged_pos_count / epi_len > 0.2:
            for i in range(3, min(len(cdr_h3) - 4, 10)):
                if random.random() < 0.4:  # 40% chance to modify
                    cdr_h3[i] = random.choice(charged_neg)
        
        # For hydrophobic epitopes, add hydrophobics and aromatics
        elif hydrophobic_count / epi_len > 0.3:
            for i in range(3, min(len(cdr_h3) - 4, 10)):
                if random.random() < 0.5:  # 50% chance to modify
                    cdr_h3[i] = random.choice(hydrophobic)
        
        template_cdrs["CDR-H3"] = ''.join(cdr_h3)
        
        # Customize other CDRs slightly
        cdr_h1 = list(template_cdrs["CDR-H1"])
        cdr_h1[4] = random.choice("NSTY")
        template_cdrs["CDR-H1"] = ''.join(cdr_h1)
        
        cdr_l3 = list(template_cdrs["CDR-L3"])
        cdr_l3[2] = random.choice("SYTN")
        template_cdrs["CDR-L3"] = ''.join(cdr_l3)
        
        return template_cdrs
    
    def _graft_cdrs_to_framework(self, cdr_sequences: Dict[str, str]) -> Tuple[str, str]:
        """Graft CDR sequences onto the framework to create full antibody chains"""
        # Get framework sequences
        heavy_framework = self.framework["heavy"]["sequence"]
        light_framework = self.framework["light"]["sequence"]
        
        # Get CDR positions in the framework
        heavy_cdrs = self.framework["heavy"].get("CDRs", {
            "CDR-H1": (26, 32),
            "CDR-H2": (52, 56),
            "CDR-H3": (95, 102)
        })
        
        light_cdrs = self.framework["light"].get("CDRs", {
            "CDR-L1": (24, 34),
            "CDR-L2": (50, 56),
            "CDR-L3": (89, 97)
        })
        
        # Insert CDRs into framework
        heavy_chain = list(heavy_framework)
        for cdr_name, (start, end) in heavy_cdrs.items():
            if cdr_name in cdr_sequences:
                # Calculate length difference
                original_len = end - start + 1
                new_len = len(cdr_sequences[cdr_name])
                
                # Replace CDR
                heavy_chain[start:end+1] = cdr_sequences[cdr_name]
                
                # Adjust other CDR positions if length changed
                if new_len != original_len:
                    diff = new_len - original_len
                    for other_cdr, (other_start, other_end) in heavy_cdrs.items():
                        if other_start > end:
                            heavy_cdrs[other_cdr] = (other_start + diff, other_end + diff)
        
        light_chain = list(light_framework)
        for cdr_name, (start, end) in light_cdrs.items():
            if cdr_name in cdr_sequences:
                # Calculate length difference
                original_len = end - start + 1
                new_len = len(cdr_sequences[cdr_name])
                
                # Replace CDR
                light_chain[start:end+1] = cdr_sequences[cdr_name]
                
                # Adjust other CDR positions if length changed
                if new_len != original_len:
                    diff = new_len - original_len
                    for other_cdr, (other_start, other_end) in light_cdrs.items():
                        if other_start > end:
                            light_cdrs[other_cdr] = (other_start + diff, other_end + diff)
        
        return ''.join(heavy_chain), ''.join(light_chain)
    
    def _extract_cdrs_from_sequences(self, heavy_chain: str, light_chain: str) -> Dict[str, str]:
        """Extract CDR sequences from full antibody sequences"""
        cdrs = {}
        
        # Get CDR positions in the framework
        heavy_cdrs = self.framework["heavy"].get("CDRs", {
            "CDR-H1": (26, 32),
            "CDR-H2": (52, 56), 
            "CDR-H3": (95, 102)
        })
        
        light_cdrs = self.framework["light"].get("CDRs", {
            "CDR-L1": (24, 34),
            "CDR-L2": (50, 56),
            "CDR-L3": (89, 97) 
        })
        
        # Extract CDRs
        for cdr_name, (start, end) in heavy_cdrs.items():
            if start < len(heavy_chain) and end < len(heavy_chain):
                cdrs[cdr_name] = heavy_chain[start:end+1]
        
        for cdr_name, (start, end) in light_cdrs.items():
            if start < len(light_chain) and end < len(light_chain):
                cdrs[cdr_name] = light_chain[start:end+1]
        
        return cdrs
    
    def _generate_cdr_variants(self, cdr_seq: str, epitope_seq: str, cdr_name: str, num_variants: int = 5) -> List[str]:
        """Generate variations of a CDR sequence"""
        variants = [cdr_seq]  # Include original
        
        # Determine positions that can be modified
        # For CDR-H3, middle positions are most variable
        if cdr_name == "CDR-H3":
            mutable_positions = list(range(2, len(cdr_seq) - 2))
        else:
            # For other CDRs, modify non-anchor residues
            mutable_positions = list(range(1, len(cdr_seq) - 1))
        
        # Generate variants
        for _ in range(num_variants - 1):
            variant = list(cdr_seq)
            
            # Make 1-3 mutations
            num_mutations = random.randint(1, min(3, len(mutable_positions)))
            positions_to_mutate = random.sample(mutable_positions, num_mutations)
            
            for pos in positions_to_mutate:
                # Determine amino acid properties to target
                # For simplicity, use basic categories
                aa_groups = {
                    "charged": "KRDEH",
                    "polar": "NQST", 
                    "hydrophobic": "AILMFWYV",
                    "small": "GASTP",
                    "aromatic": "FYW"
                }
                
                # Choose a group with bias toward complementary properties
                group_name = random.choice(list(aa_groups.keys()))
                group = aa_groups[group_name]
                
                # Choose a random aa from the group, but not the same as current
                current_aa = variant[pos]
                choices = [aa for aa in group if aa != current_aa]
                if choices:
                    variant[pos] = random.choice(choices)
            
            variants.append(''.join(variant))
        
        return variants
    
    def _score_cdr_epitope_binding(self, cdr_seq: str, epitope_seq: str) -> float:
        """Score the binding potential between a CDR and epitope"""
        # This is a simplified scoring function
        # A real implementation would use physics-based or ML-based scoring
        
        # Define complementary properties
        complementary = {
            "K": "DE",  # Lys complements Asp/Glu
            "R": "DE",  # Arg complements Asp/Glu
            "D": "KR",  # Asp complements Lys/Arg
            "E": "KR",  # Glu complements Lys/Arg
            "F": "FYWLIMV",  # Phe has hydrophobic interactions
            "Y": "FYWLIMV",  # Tyr has hydrophobic and H-bond
            "W": "FYWLIMV",  # Trp has hydrophobic interactions
            "L": "FYWLIMV",  # Leu has hydrophobic interactions
            "I": "FYWLIMV",  # Ile has hydrophobic interactions
            "M": "FYWLIMV",  # Met has hydrophobic interactions
            "V": "FYWLIMV",  # Val has hydrophobic interactions
            "S": "DENQKR",   # Ser can H-bond with charged/polar
            "T": "DENQKR",   # Thr can H-bond with charged/polar
            "N": "DENQKRST", # Asn can H-bond with charged/polar
            "Q": "DENQKRST"  # Gln can H-bond with charged/polar
        }
        
        # Count complementary pairs
        complement_count = 0
        for cdr_aa in cdr_seq:
            if cdr_aa in complementary:
                for epitope_aa in epitope_seq:
                    if epitope_aa in complementary[cdr_aa]:
                        complement_count += 1
        
        # Normalize by length
        score = complement_count / (len(cdr_seq) * len(epitope_seq))
        
        # Add some randomness to simulate uncertainty in real scoring
        score = score * (0.9 + 0.2 * random.random())
        
        return score
    
    def _predict_binding_energy(self, heavy_chain: str, light_chain: str, epitope_seq: str) -> float:
        """Predict binding energy between antibody and epitope"""
        # This is a simplified model
        # In a real implementation, this would use molecular modeling or ML
        
        # Extract CDRs
        cdrs = self._extract_cdrs_from_sequences(heavy_chain, light_chain)
        
        # Score each CDR
        scores = {}
        weights = {
            "CDR-H1": 0.15,
            "CDR-H2": 0.20,
            "CDR-H3": 0.40,  # H3 is most important
            "CDR-L1": 0.10,
            "CDR-L2": 0.05,
            "CDR-L3": 0.10
        }
        
        for cdr_name, cdr_seq in cdrs.items():
            scores[cdr_name] = self._score_cdr_epitope_binding(cdr_seq, epitope_seq)
        
        # Weighted sum of scores
        binding_score = sum(weights.get(cdr, 0.1) * scores.get(cdr, 0) for cdr in cdrs)
        
        # Convert to energy (negative values = stronger binding)
        # Scale to realistic range of binding energies
        binding_energy = -10.0 * binding_score - 5.0 * random.random()
        
        return round(binding_energy, 2)
    
    def _convert_energy_to_affinity(self, energy: float) -> float:
        """Convert binding energy to affinity in nanomolar"""
        # Using a simplified relationship between ΔG and Kd
        # ΔG = -RT ln(Kd) => Kd = exp(-ΔG/RT)
        # For room temperature, RT ≈ 0.6 kcal/mol
        
        # Convert energy to affinity in nM
        # Lower (more negative) energy = stronger binding = lower Kd
        kd_nm = 1000.0 * np.exp(energy / 0.6)
        
        # Ensure reasonable affinity range (1 pM to 10 μM)
        kd_nm = max(0.001, min(10000.0, kd_nm))
        
        return round(kd_nm, 2)
    
    def _calculate_developability(self, heavy_chain: str, light_chain: str) -> Dict[str, float]:
        """Calculate developability metrics for an antibody"""
        # This is a simplified implementation
        # Real implementations would use more sophisticated models
        
        # Check for developability issues
        issues = []
        
        # 1. Check for unusual amino acid frequencies
        aa_freq = {}
        for aa in heavy_chain + light_chain:
            aa_freq[aa] = aa_freq.get(aa, 0) + 1
        total_aa = len(heavy_chain) + len(light_chain)
        
        # Check for high cysteine content (unpaired cysteines can cause problems)
        cys_freq = aa_freq.get('C', 0) / total_aa
        if cys_freq > 0.05:  # More than 5% Cys
            issues.append("high_cysteine")
        
        # 2. Check for hydrophobicity of CDRs
        hydrophobic = "AVILMFYW"
        cdrs = self._extract_cdrs_from_sequences(heavy_chain, light_chain)
        
        cdr_h3 = cdrs.get("CDR-H3", "")
        if cdr_h3:
            h3_hydrophobicity = sum(1 for aa in cdr_h3 if aa in hydrophobic) / len(cdr_h3)
            if h3_hydrophobicity > 0.6:  # More than 60% hydrophobic
                issues.append("high_hydrophobicity_cdr_h3")
        
        # 3. Assign a developability score
        base_score = 0.90  # Start with a good score
        deductions = {
            "high_cysteine": 0.15,
            "high_hydrophobicity_cdr_h3": 0.10
        }
        
        quality_score = base_score - sum(deductions.get(issue, 0) for issue in issues)
        quality_score = max(0.0, min(1.0, quality_score))  # Ensure 0-1 range
        
        # Manufacturability adds some variability
        manufacturability = 0.85 + 0.15 * random.random()
        
        # Overall developability
        developability = (quality_score * 0.7 + manufacturability * 0.3)
        
        return {
            "overall_score": round(developability, 2),
            "quality_score": round(quality_score, 2),
            "manufacturability": round(manufacturability, 2),
            "issues": issues
        }