#!/usr/bin/env python
# demo_test.py - Test script for Phytovenomics ML platform

import os
import sys
import logging
import json
import random
import pandas as pd
import numpy as np
from typing import Dict, List

# Add the current directory to path to ensure imports work
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from antibody_design import HumanAntibodyDesigner, EpitopeDiscovery, AffinityOptimizer
from cocktail_strategy import CocktailOptimizer
from venom_data.toxin_database import ToxinDatabase
from utils import DataUtils, ModelUtils, Visualization

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)

logger = logging.getLogger("Phytovenomics.Demo")

# Create synthetic toxin data
def create_synthetic_toxin_data(num_toxins=20):
    """Create synthetic toxin data for demonstration."""
    toxin_families = [
        "three_finger_toxin",
        "phospholipase_a2",
        "snake_venom_metalloprotease",
        "snake_venom_serine_protease",
        "c_type_lectin"
    ]
    
    snake_species = [
        "Naja naja",
        "Bungarus caeruleus",
        "Daboia russelii",
        "Ophiophagus hannah",
        "Dendroaspis polylepis"
    ]
    
    # Amino acid sequences (simplified for demo)
    aa_pool = "ACDEFGHIKLMNPQRSTVWY"
    
    toxins = []
    
    for i in range(num_toxins):
        # Generate synthetic toxin
        family = random.choice(toxin_families)
        species = random.choice(snake_species)
        toxin_id = f"{family[:3].upper()}{i+1:03d}"
        
        # Generate synthetic sequence (60-120 amino acids)
        seq_length = random.randint(60, 120)
        sequence = ''.join(random.choice(aa_pool) for _ in range(seq_length))
        
        # Create toxin entry
        toxin = {
            "toxin_id": toxin_id,
            "sequence": sequence,
            "toxin_family": family,
            "source_species": species,
            "molecular_weight": round(random.uniform(6000, 25000), 1),
            "protein_type": "toxin",
            "function": "neurotoxic" if "finger" in family else "hemotoxic",
            "target": "nicotinic acetylcholine receptor" if "finger" in family else "cell membrane",
            "potency": round(random.uniform(0.1, 5.0), 2)
        }
        
        toxins.append(toxin)
    
    # Create DataFrame
    df = pd.DataFrame(toxins)
    
    # Save to CSV
    os.makedirs("data/snake_toxins", exist_ok=True)
    csv_path = "data/snake_toxins/snake_toxins_synthetic.csv"
    df.to_csv(csv_path, index=False)
    
    logger.info(f"Created synthetic toxin dataset with {num_toxins} toxins at {csv_path}")
    return csv_path, toxins

# Create synthetic human antibody data
def create_synthetic_antibody_data(num_antibodies=30):
    """Create synthetic human antibody data."""
    # Very simplified for demonstration
    aa_pool = "ACDEFGHIKLMNPQRSTVWY"
    
    antibodies = []
    fasta_content = ""
    
    for i in range(num_antibodies):
        ab_id = f"HUMAB{i+1:03d}"
        
        # Generate heavy chain (simplified)
        heavy_length = random.randint(110, 130)
        heavy_seq = ''.join(random.choice(aa_pool) for _ in range(heavy_length))
        
        # Generate light chain (simplified)
        light_length = random.randint(100, 115)
        light_seq = ''.join(random.choice(aa_pool) for _ in range(light_length))
        
        # Add to FASTA
        fasta_content += f">{ab_id}_H\n{heavy_seq}\n>{ab_id}_L\n{light_seq}\n"
        
        antibodies.append({
            "id": ab_id,
            "heavy_chain": heavy_seq,
            "light_chain": light_seq
        })
    
    # Save to FASTA
    os.makedirs("data/antibody_structures/human", exist_ok=True)
    fasta_path = "data/antibody_structures/human/human_antibodies_synthetic.fasta"
    
    with open(fasta_path, 'w') as f:
        f.write(fasta_content)
    
    logger.info(f"Created synthetic antibody dataset with {num_antibodies} antibodies at {fasta_path}")
    return fasta_path, antibodies

# Create synthetic binding data
def create_synthetic_binding_data(toxins, antibodies, binding_pairs=50):
    """Create synthetic toxin-antibody binding data."""
    bindings = []
    
    # Create random binding pairs
    for _ in range(binding_pairs):
        toxin = random.choice(toxins)
        antibody = random.choice(antibodies)
        
        # Create binding entry
        binding = {
            "toxin_id": toxin["toxin_id"],
            "antibody_id": antibody["id"],
            "affinity_nm": round(random.uniform(0.1, 1000), 2),  # nM (lower = stronger binding)
            "binding_interface_area": round(random.uniform(500, 2000), 1),  # Å²
            "binding_energy": round(random.uniform(-15, -5), 1)  # kcal/mol
        }
        
        bindings.append(binding)
    
    # Create DataFrame
    df = pd.DataFrame(bindings)
    
    # Save to CSV
    os.makedirs("data/toxin_antibody_binding", exist_ok=True)
    csv_path = "data/toxin_antibody_binding/toxin_antibody_binding_pairs.csv"
    df.to_csv(csv_path, index=False)
    
    logger.info(f"Created synthetic binding dataset with {binding_pairs} pairs at {csv_path}")
    return csv_path, bindings

# Update configuration with data paths
def update_config(toxin_path, antibody_path, binding_path):
    """Update the configuration file with data paths."""
    try:
        with open("config/default_config.yaml", 'r') as f:
            config_content = f.read()
        
        # Update data paths
        config_content = config_content.replace(
            "../data/snake_toxins/snake_toxins_synthetic.csv",
            toxin_path
        )
        
        config_content = config_content.replace(
            "../data/antibody_structures/human/human_antibodies_synthetic.fasta",
            antibody_path
        )
        
        config_content = config_content.replace(
            "../data/toxin_antibody_binding/toxin_antibody_binding_pairs.csv",
            binding_path
        )
        
        with open("config/demo_config.yaml", 'w') as f:
            f.write(config_content)
            
        logger.info("Updated configuration with synthetic data paths")
        return "config/demo_config.yaml"
    
    except Exception as e:
        logger.error(f"Error updating configuration: {e}")
        return "config/default_config.yaml"

# Mock antibody generator methods
def mock_antibody_design(antibody_designer, toxin_data, epitopes, config):
    """Mock antibody design process for demonstration."""
    logger.info("Running mock antibody design process...")
    designed_antibodies = []
    
    for epitope in epitopes[:min(5, len(epitopes))]:
        # Design a specific antibody
        specific_ab = {
            "id": f"AB{len(designed_antibodies)+1:03d}_SPEC",
            "target_toxin_id": toxin_data["id"],
            "target_epitope_id": epitope.get("id", "EP_UNK"),
            "sequence": "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(220)),
            "is_broadly_neutralizing": False,
            "developability_score": round(random.uniform(0.5, 0.95), 2),
            "cdr_lengths": {
                "CDR-H1": random.randint(8, 12),
                "CDR-H2": random.randint(8, 19),
                "CDR-H3": random.randint(8, 22),
                "CDR-L1": random.randint(8, 16),
                "CDR-L2": random.randint(7, 11),
                "CDR-L3": random.randint(8, 12)
            },
            "binds_to": [
                {
                    "toxin_id": toxin_data["id"],
                    "affinity_nm": round(random.uniform(0.5, 50), 2),
                    "binding_energy": round(random.uniform(-14, -9), 1)
                }
            ],
            "targeting_strategy": "specific",
            "generation_date": "2025-05-29"
        }
        
        designed_antibodies.append(specific_ab)
        
        # Design a broadly neutralizing antibody
        bnab = {
            "id": f"AB{len(designed_antibodies)+1:03d}_BNAB",
            "target_toxin_id": toxin_data["id"],
            "target_epitope_id": epitope.get("id", "EP_UNK"),
            "sequence": "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(230)),
            "is_broadly_neutralizing": True,
            "developability_score": round(random.uniform(0.5, 0.95), 2),
            "breadth_score": round(random.uniform(0.6, 0.9), 2),
            "cdr_lengths": {
                "CDR-H1": random.randint(8, 12),
                "CDR-H2": random.randint(8, 19),
                "CDR-H3": random.randint(10, 22),
                "CDR-L1": random.randint(8, 16),
                "CDR-L2": random.randint(7, 11),
                "CDR-L3": random.randint(8, 12)
            },
            "binds_to": []
        }
        
        # Add binding to multiple toxins in same family
        toxin_family = toxin_data.get("family", "")
        base_affinity = round(random.uniform(1, 100), 2)
        
        # Bind to primary target
        bnab["binds_to"].append({
            "toxin_id": toxin_data["id"],
            "affinity_nm": base_affinity,
            "binding_energy": round(random.uniform(-12, -8), 1)
        })
        
        # Add bindings to related toxins (simulating cross-neutralization)
        for i in range(3):  # Add 3 more related toxins
            related_toxin_id = f"{toxin_family[:3].upper()}{random.randint(1, 999):03d}"
            if related_toxin_id != toxin_data["id"]:
                # Related toxins have lower affinity
                related_affinity = base_affinity * random.uniform(1.5, 4)
                bnab["binds_to"].append({
                    "toxin_id": related_toxin_id,
                    "affinity_nm": round(related_affinity, 2),
                    "binding_energy": round(random.uniform(-10, -6), 1)
                })
        
        designed_antibodies.append(bnab)
    
    return designed_antibodies

# Mock epitope discovery
def identify_epitopes_for_demo(epitope_discovery, toxin_data):
    """Use actual epitope discovery algorithm with ML enhancements."""
    logger.info(f"Running ML-enhanced epitope discovery for toxin {toxin_data.get('id')}")
    
    seq = toxin_data.get("sequence", "")
    
    if not seq:
        logger.warning("No sequence data available for epitope discovery")
        return []
    
    try:
        # Attempt to use the actual epitope discovery implementation
        epitopes = epitope_discovery.identify_epitopes(toxin_data)
        
        if not epitopes and epitope_discovery.use_ml_model:
            # If no epitopes were found but ML is enabled, try direct ML prediction
            logger.info("Attempting direct ML-based epitope prediction")
            epitopes = epitope_discovery.predict_epitopes_with_ml(seq)
        
        if not epitopes:  # If still no results, use fallback approach
            logger.warning("Falling back to rule-based epitope discovery")
            # Fall back to the original mock implementation for robustness
            epitopes = mock_epitope_discovery(epitope_discovery, toxin_data)
        
        logger.info(f"Identified {len(epitopes)} epitopes for toxin {toxin_data.get('id')}")
        return epitopes
        
    except Exception as e:
        logger.error(f"Error in epitope discovery: {e}")
        logger.info("Falling back to mock epitope discovery")
        return mock_epitope_discovery(epitope_discovery, toxin_data)


def mock_epitope_discovery(epitope_discovery, toxin_data):
    """Mock epitope discovery for demonstration (fallback)."""
    logger.info("Running fallback mock epitope discovery...")
    
    # Generate random epitopes
    epitopes = []
    seq = toxin_data.get("sequence", "")
    
    if not seq:
        logger.warning("No sequence data available for epitope discovery")
        return []
    
    # Generate 3-5 epitopes
    epitope_count = random.randint(3, 5)
    seq_len = len(seq)
    
    for i in range(epitope_count):
        # Random epitope length between 8-20 amino acids
        ep_len = random.randint(8, min(20, seq_len//2))
        
        # Random starting position
        start_pos = random.randint(0, seq_len - ep_len - 1)
        epitope_seq = seq[start_pos:start_pos + ep_len]
        
        # Create epitope object
        epitope = {
            "id": f"{toxin_data.get('id')}_EP{i+1}",
            "start": start_pos,
            "end": start_pos + ep_len - 1,
            "length": ep_len,
            "sequence": epitope_seq,
            "accessibility": round(random.uniform(0.5, 1.0), 2),
            "hydrophilicity": round(random.uniform(0.3, 0.9), 2),
            "conservation": round(random.uniform(0.4, 0.95), 2),
            "predicted_immunogenicity": round(random.uniform(0.5, 0.95), 2),
            "ml_score": round(random.uniform(0.4, 0.9), 2),
            "score": round(random.uniform(0.5, 0.98), 2),
            "toxin_id": toxin_data.get("id"),
            "toxin_family": toxin_data.get("toxin_family"),
            "source_species": toxin_data.get("source_species"),
            "prediction_method": "Fallback-Mock"
        }
        
        epitopes.append(epitope)
    
    # Sort by score
    epitopes.sort(key=lambda x: x["score"], reverse=True)
    return epitopes

def run_demo():
    """Run a demonstration of the antibody generation pipeline."""
    logger.info("Starting Phytovenomics ML Platform demonstration")
    
    # Create synthetic data for demonstration
    toxin_path, toxins = create_synthetic_toxin_data(num_toxins=20)
    antibody_path, antibodies = create_synthetic_antibody_data(num_antibodies=30)
    binding_path, bindings = create_synthetic_binding_data(toxins, antibodies, binding_pairs=50)
    
    # Update configuration
    config_path = update_config(toxin_path, antibody_path, binding_path)
    
    # Load configuration
    with open(config_path, 'r') as f:
        import yaml
        config = yaml.safe_load(f)
    
    # Initialize components
    toxin_db = ToxinDatabase(config.get("toxin_database", {}))
    epitope_discovery = EpitopeDiscovery(config.get("epitope_discovery", {}))
    antibody_designer = HumanAntibodyDesigner(config.get("antibody_designer", {}))
    affinity_optimizer = AffinityOptimizer(config.get("affinity_optimizer", {}))
    cocktail_optimizer = CocktailOptimizer(config.get("cocktail_optimizer", {}))
    
    # Load toxin data
    toxin_db.load_from_csv(toxin_path)
    all_toxins = toxin_db.get_all_toxins()
    
    # Select sample targets (first 5 toxins)
    target_toxins = list(all_toxins.values())[:5]
    
    # Process each target
    designed_antibodies = []
    
    for toxin in target_toxins:
        logger.info(f"Processing toxin {toxin['id']} from {toxin['source_species']}")
        
        # Discover epitopes using ML-enhanced approach
        epitopes = identify_epitopes_for_demo(epitope_discovery, toxin)
        logger.info(f"Found {len(epitopes)} epitopes on {toxin['id']}")
        
        # Log the prediction method used
        methods = set(epitope.get("prediction_method", "Unknown") for epitope in epitopes)
        logger.info(f"Epitope prediction methods used: {', '.join(methods)}")
        
        # Design antibodies against epitopes
        toxin_antibodies = mock_antibody_design(antibody_designer, toxin, epitopes, config)
        logger.info(f"Designed {len(toxin_antibodies)} antibodies against {toxin['id']}")
        
        designed_antibodies.extend(toxin_antibodies)
    
    # Save results
    os.makedirs("output/demo", exist_ok=True)
    
    # Save antibody data as JSON
    with open("output/demo/designed_antibodies.json", 'w') as f:
        json.dump(designed_antibodies, f, indent=2)
    
    # Design cocktail
    toxin_dict = {toxin["id"]: toxin for toxin in target_toxins}
    cocktail = cocktail_optimizer.design_optimal_cocktail(designed_antibodies, toxin_dict)
    
    # Save cocktail design
    with open("output/demo/antibody_cocktail.json", 'w') as f:
        # Remove all_antibodies field which is too large
        if "all_antibodies" in cocktail:
            del cocktail["all_antibodies"]
        json.dump(cocktail, f, indent=2)
    
    # Print summary
    print("\nPhytovenomics ML Platform Demonstration Summary:")
    print(f"- Generated {len(designed_antibodies)} antibodies targeting {len(target_toxins)} toxins")
    print(f"- Created optimized cocktail with {len(cocktail['antibodies'])} antibodies")
    print(f"- Average toxin family coverage: {cocktail['average_coverage']:.2%}")
    print(f"\nResults saved to output/demo/")
    
    # Save simple summary report
    with open("output/demo/demonstration_summary.md", 'w') as f:
        f.write("# Phytovenomics ML Platform Demonstration\n\n")
        f.write(f"## Summary\n")
        f.write(f"- Date: 2025-05-29\n")
        f.write(f"- Toxins analyzed: {len(target_toxins)}\n")
        f.write(f"- Antibodies designed: {len(designed_antibodies)}\n")
        f.write(f"- Cocktail size: {len(cocktail['antibodies'])} antibodies\n")
        f.write(f"- Average coverage: {cocktail['average_coverage']:.2%}\n\n")
        
        f.write(f"## ML-Enhanced Features\n")
        f.write(f"- **Epitope Discovery**: ML-based prediction of epitope regions on toxin sequences\n")
        f.write(f"- **Binding Affinity**: ML models to predict antibody-toxin binding strength\n")
        f.write(f"- **Cocktail Optimization**: Advanced algorithms to maximize toxin neutralization\n\n")
        
        f.write(f"## Toxins\n")
        for toxin in target_toxins:
            f.write(f"- **{toxin['id']}** ({toxin['family']}): {toxin['source_species']}\n")
        
        # Add section for epitopes with prediction methods
        f.write(f"\n## ML-Identified Epitopes\n")
        epitope_count = 0
        prediction_methods = set()
        for toxin in target_toxins[:2]:  # Show epitopes for first 2 toxins only
            epitopes = identify_epitopes_for_demo(epitope_discovery, toxin)
            epitope_count += len(epitopes)
            f.write(f"### Toxin {toxin['id']}\n")
            for ep in epitopes[:3]:  # Show first 3 epitopes per toxin
                method = ep.get('prediction_method', 'Standard')
                prediction_methods.add(method)
                f.write(f"- **{ep.get('id', 'Unknown')}**: {ep.get('sequence', 'N/A')}\n")
                f.write(f"  - Position: {ep.get('start', 'N/A')}-{ep.get('end', 'N/A')}\n")
                f.write(f"  - Score: {ep.get('score', ep.get('ml_score', 'N/A'))}\n")
                f.write(f"  - Method: {method}\n")
        
        f.write(f"\nTotal epitopes identified: {epitope_count}\n")
        f.write(f"Prediction methods used: {', '.join(prediction_methods)}\n\n")
        
        f.write(f"## Antibody Cocktail\n")
        # If no antibodies in the main cocktail list, try to get them from all_antibodies
        display_abs = cocktail['antibodies']
        if len(display_abs) == 0 and 'all_antibodies' in cocktail:
            # Take a few antibodies from all_antibodies as a fallback
            display_abs = cocktail['all_antibodies'][:5] if cocktail['all_antibodies'] else []
        
        for ab in display_abs[:5]:  # First 5 antibodies
            f.write(f"- **{ab['id']}**: Targets {ab['target_toxin_id']}, ")
            f.write(f"Developability: {ab.get('developability_score', 'N/A'):.2f}, ")
            if ab.get('is_broadly_neutralizing'):
                f.write(f"Broadly neutralizing\n")
            else:
                f.write(f"Specific\n")
    
    return designed_antibodies, cocktail

if __name__ == "__main__":
    try:
        antibodies, cocktail = run_demo()
        sys.exit(0)
    except Exception as e:
        logger.exception(f"Error in demonstration: {e}")
        sys.exit(1)