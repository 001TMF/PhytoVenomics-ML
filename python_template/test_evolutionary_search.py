#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test script for the evolutionary search module with RosettaFold integration.
"""

import os
import sys
import logging
import argparse
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Dict, Optional, Tuple

# Add parent directory to path to import modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from antibody_design.evolutionary_search import EvolutionarySearch, Individual, RosettaFoldClient
from antibody_design.cdr_processor import CDRProcessor
from src.structure_prediction.hybrid_prediction import HybridPredictionPipeline
from utils.visualization import plot_antibody_evolution, plot_antibody_structures

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/evolutionary_search_test.log"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger("evolutionary_search_test")


def load_test_antibodies(num_samples: int = 5) -> List[Individual]:
    """
    Load test antibody sequences and create Individual instances.
    
    Args:
        num_samples: Number of antibody samples to load
        
    Returns:
        List of Individual instances
    """
    # Sample antibody sequences (these are real antibody sequences)
    antibodies = [
        {
            "name": "Trastuzumab",
            "heavy": "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVRGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS",
            "light": "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
            "cdrs": {
                "H1": "GFNIKDTYIH",
                "H2": "ARIYPTNGYTRYADS",
                "H3": "SRWGGDGFYAMDY",
                "L1": "RASQDVNTAVA",
                "L2": "SASFLY",
                "L3": "QQHYTTPPT"
            }
        },
        {
            "name": "Adalimumab",
            "heavy": "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS",
            "light": "DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIK",
            "cdrs": {
                "H1": "GFTFDDYAMH",
                "H2": "AITWNSGHIDYADS",
                "H3": "AKVSYLSTASSLDY",
                "L1": "RASQGIRNY",
                "L2": "AASTLQS",
                "L3": "QRYNRAPYT"
            }
        },
        {
            "name": "Pembrolizumab",
            "heavy": "QVQLVQSGVEVKKPGASVKVSCKASGYTFTDYYMYWVRQAPGQGLEWMGDINPSNGGTNYNEKFKGRVTLTTDTSTSTAYMELRSLRSDDTAVYYCARWGGDGFYAMDYWGQGTLVTVSS",
            "light": "EIVLTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYDASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSLTFGGGTKVEIK",
            "cdrs": {
                "H1": "GYTFTDYYMY",
                "H2": "DINPSNGGTNYNEK",
                "H3": "ARWGGDGFYAMDY",
                "L1": "RASQSVSSYLA",
                "L2": "DASSRAT",
                "L3": "QQYGSSLT"
            }
        },
        {
            "name": "Bevacizumab",
            "heavy": "EVQLVESGGGLVQPGGSLRLSCAASGYTFTNYGMNWVRQAPGKGLEWVGWINTYTGEPTYAADFKRRFTFSLDTSKSTAYLQMNSLRAEDTAVYYCAKYPHYYGSSHWYFDVWGQGTLVTVSS",
            "light": "DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYLYAPYTFGQGTKVEIK",
            "cdrs": {
                "H1": "GYTFTNYGMN",
                "H2": "WINTYTGEPTYAAD",
                "H3": "AKYPHYYGSSHWYFDY",
                "L1": "RASQDVSTAVA",
                "L2": "SASFLYSG",
                "L3": "QQYLYAPYT"
            }
        },
        {
            "name": "Nivolumab",
            "heavy": "QVQLVESGGGVVQPGRSLRLSCAASGFTFSRQMHWVRQAPGKGLEWVAVITYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAREGNYYGSGSPYFDYWGQGTLVTVSS",
            "light": "EIVLTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQRSNWPPTFGQGTKVEIK",
            "cdrs": {
                "H1": "GFTFSRQMH",
                "H2": "VITYDGSNKYYA",
                "H3": "AREGNYYGSGSPYFDY",
                "L1": "RASQSVSSYLA",
                "L2": "DASNRAT",
                "L3": "QQRSNWPPT"
            }
        }
    ]
    
    # Create Individual instances
    individuals = []
    for i in range(min(num_samples, len(antibodies))):
        ab = antibodies[i]
        individuals.append(Individual(
            heavy_chain=ab["heavy"],
            light_chain=ab["light"],
            cdr_sequences=ab["cdrs"]
        ))
        # Set the name as target_id since name is not directly accepted
        individuals[-1].target_id = ab["name"]
    
    return individuals


def load_test_epitopes(num_samples: int = 1) -> List[Dict]:
    """
    Load test epitope sequences and related data.
    
    Args:
        num_samples: Number of epitope samples to load
        
    Returns:
        List of epitope dictionaries
    """
    epitopes = [
        {
            "id": "spike_rbd_001", 
            "sequence": "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLD",
            "structure": None,
            "organism": "SARS-CoV-2",
            "region": "RBD"
        },
        {
            "id": "snake_toxin_001",
            "sequence": "LKCNQLIPPFWKTCPEGKNLCYKMMLASKKMVPVKRGCINVCPKNSALVKYVCCNTDRCN",
            "structure": None,
            "organism": "Dendroaspis polylepis",
            "region": "Neurotoxin"
        },
    ]
    
    return epitopes[:num_samples]


def run_standard_optimization_test(config_path: str) -> Tuple[List[Individual], Dict]:
    """
    Run a standard optimization test using the evolutionary search algorithm.
    
    Args:
        config_path: Path to the configuration file
        
    Returns:
        Tuple of optimized antibodies and test results
    """
    logger.info("Starting standard optimization test")
    
    # Import config loader
    from antibody_design.config_loader import load_config
    
    # Load configuration and initialize evolutionary search
    config = load_config(config_path)
    search = EvolutionarySearch(config)
    
    # Load test data
    initial_population = load_test_antibodies(5)
    target_epitope = load_test_epitopes(1)[0]
    
    logger.info(f"Initial population: {len(initial_population)} antibodies")
    logger.info(f"Target epitope: {target_epitope['id']} ({target_epitope['organism']})")
    
    # Run optimization
    results = search.optimize_antibody(initial_population, target_epitope)
    
    # Print results
    logger.info(f"Optimization complete! Found {len(results)} optimal antibodies.")
    
    for i, antibody in enumerate(results):
        logger.info(f"\nOptimal Antibody {i+1}:")
        logger.info(f"ID: {antibody.target_id if hasattr(antibody, 'target_id') else 'Unknown'}")
        logger.info(f"Binding Energy: {antibody.binding_energy:.2f}")
        logger.info(f"Stability: {antibody.stability:.2f}")
        logger.info(f"Developability: {antibody.developability:.2f}")
        logger.info(f"Manufacturability: {antibody.manufacturability:.2f}")
        
        # Print CDR sequences
        logger.info("CDR sequences:")
        for cdr_name, cdr_seq in antibody.cdr_sequences.items():
            logger.info(f"  {cdr_name}: {cdr_seq}")
    
    # Prepare test results
    test_results = {
        "num_initial": len(initial_population),
        "num_optimal": len(results),
        "target_epitope": target_epitope["id"],
        "best_binding": min([ab.binding_energy for ab in results]) if results else None,
        "best_stability": max([ab.stability for ab in results]) if results else None,
        "best_overall": results[0].name if results else None,
    }
    
    return results, test_results


def run_hybrid_structure_test(config_path: str, optimized_antibodies: List[Individual]) -> Dict:
    """
    Test integration with the hybrid structure prediction pipeline.
    
    Args:
        config_path: Path to the configuration file
        optimized_antibodies: List of optimized antibodies from the evolutionary search
        
    Returns:
        Dictionary of test results
    """
    logger.info("Starting hybrid structure prediction test")
    
    try:
        # Initialize hybrid prediction pipeline
        hybrid_pipeline = HybridPredictionPipeline(
            use_mock=(os.getenv("USE_MOCK_PREDICTION", "True").lower() == "true")
        )
        
        # Test structure prediction for the top antibody
        if optimized_antibodies:
            top_antibody = optimized_antibodies[0]
            
            logger.info(f"Predicting structure for top antibody: {top_antibody.name}")
            
            # Combine heavy and light chains
            fasta_content = f">Top_AB_H\n{top_antibody.heavy_chain}\n>Top_AB_L\n{top_antibody.light_chain}"
            
            # Create temp fasta file
            output_dir = Path("output/evolutionary_search_test")
            os.makedirs(output_dir, exist_ok=True)
            
            fasta_path = output_dir / "top_antibody.fasta"
            with open(fasta_path, 'w') as f:
                f.write(fasta_content)
            
            # Predict structure
            result = hybrid_pipeline.predict_structure(
                fasta_path=str(fasta_path),
                output_dir=str(output_dir)
            )
            
            # Log results
            if result.get("success", False):
                logger.info(f"Structure prediction successful.")
                logger.info(f"PDB file: {result.get('pdb_path')}")
                logger.info(f"Prediction method: {result.get('method_used')}")
                logger.info(f"Confidence score: {result.get('confidence', 'N/A')}")
            else:
                logger.error(f"Structure prediction failed: {result.get('error', 'Unknown error')}")
                
            return {
                "success": result.get("success", False),
                "method_used": result.get("method_used"),
                "confidence": result.get("confidence"),
                "pdb_path": result.get("pdb_path"),
            }
        else:
            logger.warning("No optimized antibodies available for structure prediction test")
            return {"success": False, "error": "No optimized antibodies available"}
            
    except Exception as e:
        logger.exception(f"Error during hybrid structure prediction: {str(e)}")
        return {"success": False, "error": str(e)}


def run_multi_objective_test(config_path: str) -> Dict:
    """
    Test multi-objective optimization with various weight combinations.
    
    Args:
        config_path: Path to the configuration file
        
    Returns:
        Dictionary of test results
    """
    logger.info("Starting multi-objective optimization test")
    
    # Weight configurations to test
    weight_configs = [
        {"name": "Balanced", "weights": {"binding": 1.0, "stability": 1.0, "developability": 1.0, "manufacturability": 1.0}},
        {"name": "Binding Focus", "weights": {"binding": 2.0, "stability": 0.5, "developability": 0.5, "manufacturability": 0.5}},
        {"name": "Stability Focus", "weights": {"binding": 0.5, "stability": 2.0, "developability": 0.5, "manufacturability": 0.5}},
        {"name": "Developability Focus", "weights": {"binding": 0.5, "stability": 0.5, "developability": 2.0, "manufacturability": 0.5}},
    ]
    
    results = {}
    
    for weight_config in weight_configs:
        logger.info(f"Testing weight configuration: {weight_config['name']}")
        
        # Create new config with updated weights
        search = EvolutionarySearch(config_path)
        search.objective_weights = weight_config["weights"]
        
        # Load test data
        initial_population = load_test_antibodies(3)  # Use fewer samples for faster tests
        target_epitope = load_test_epitopes(1)[0]
        
        # Run optimization with fewer generations
        search.config["num_generations"] = 5  # Reduced for testing
        optimized = search.optimize_antibody(initial_population, target_epitope)
        
        if optimized:
            top_antibody = optimized[0]
            results[weight_config["name"]] = {
                "binding_energy": top_antibody.binding_energy,
                "stability": top_antibody.stability,
                "developability": top_antibody.developability,
                "manufacturability": top_antibody.manufacturability,
                "name": top_antibody.name
            }
        else:
            results[weight_config["name"]] = {"error": "No results found"}
    
    # Log comparison
    logger.info("\nMulti-objective optimization results:")
    for config_name, res in results.items():
        logger.info(f"\n{config_name}:")
        for key, value in res.items():
            if isinstance(value, float):
                logger.info(f"  {key}: {value:.2f}")
            else:
                logger.info(f"  {key}: {value}")
    
    return results


def main():
    """
    Main function to run the evolutionary search tests.
    """
    parser = argparse.ArgumentParser(description='Test Evolutionary Search Algorithm')
    parser.add_argument('--config', default='config/evolutionary_search_config.yaml',
                       help='Path to configuration file')
    parser.add_argument('--test-type', choices=['standard', 'hybrid', 'multi', 'all'],
                       default='all', help='Type of test to run')
    parser.add_argument('--output-dir', default='output/evolutionary_search_test',
                       help='Directory to save test outputs')
    parser.add_argument('--mock', action='store_true', 
                       help='Use mock implementations instead of actual models')
    
    args = parser.parse_args()
    
    # Set mock mode environment variable
    os.environ["USE_MOCK_PREDICTION"] = "True" if args.mock else "False"
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    config_path = args.config
    logger.info(f"Using configuration file: {config_path}")
    
    # Run selected tests
    if args.test_type in ['standard', 'all']:
        optimized_antibodies, standard_results = run_standard_optimization_test(config_path)
        
        # Save results
        with open(os.path.join(args.output_dir, "standard_test_results.txt"), 'w') as f:
            for key, value in standard_results.items():
                f.write(f"{key}: {value}\n")
    else:
        optimized_antibodies = []
    
    if args.test_type in ['hybrid', 'all'] and optimized_antibodies:
        hybrid_results = run_hybrid_structure_test(config_path, optimized_antibodies)
        
        # Save results
        with open(os.path.join(args.output_dir, "hybrid_test_results.txt"), 'w') as f:
            for key, value in hybrid_results.items():
                f.write(f"{key}: {value}\n")
    
    if args.test_type in ['multi', 'all']:
        multi_results = run_multi_objective_test(config_path)
        
        # Save results
        with open(os.path.join(args.output_dir, "multi_objective_test_results.txt"), 'w') as f:
            f.write("Multi-objective optimization results:\n")
            for config_name, res in multi_results.items():
                f.write(f"\n{config_name}:\n")
                for key, value in res.items():
                    f.write(f"  {key}: {value}\n")
    
    logger.info("All tests completed!")


if __name__ == "__main__":
    main()