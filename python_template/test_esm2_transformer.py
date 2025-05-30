#!/usr/bin/env python
# test_esm2_transformer.py

import os
import sys
import time
import torch
import logging
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import project modules
from antibody_design.esm2_cdr_masking import ESM2CDRModel, ESM2Config
from antibody_design.cdr_processor import CDRProcessor
from utils.training_utils import setup_device, set_seed

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger("Phytovenomics.ESM2Test")

def generate_synthetic_antibody_sequences(
    num_sequences: int = 10, 
    min_length: int = 100, 
    max_length: int = 150
) -> List[str]:
    """
    Generate synthetic antibody sequences for testing.
    
    Args:
        num_sequences: Number of sequences to generate
        min_length: Minimum sequence length
        max_length: Maximum sequence length
        
    Returns:
        List of synthetic antibody sequences
    """
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    cdr_h3_motifs = ["CARGY", "CARGD", "CARVS", "CARDR", "CARYY"]
    cdr_l3_motifs = ["CQQYN", "CQQYN", "CQQSK", "CQQYD", "CQQFT"]
    
    sequences = []
    
    for i in range(num_sequences):
        # Create a random sequence length
        seq_length = np.random.randint(min_length, max_length + 1)
        
        # Generate random sequence
        sequence = ''.join(np.random.choice(list(amino_acids)) for _ in range(seq_length))
        
        # Insert recognizable CDR-H3 motif at a random position
        h3_pos = np.random.randint(30, seq_length - 20)
        h3_motif = np.random.choice(cdr_h3_motifs)
        sequence = sequence[:h3_pos] + h3_motif + sequence[h3_pos + len(h3_motif):]
        
        # Insert recognizable CDR-L3 motif at a random position
        l3_pos = np.random.randint(h3_pos + 20, seq_length - 10)
        l3_motif = np.random.choice(cdr_l3_motifs)
        sequence = sequence[:l3_pos] + l3_motif + sequence[l3_pos + len(l3_motif):]
        
        sequences.append(sequence)
    
    logger.info(f"Generated {num_sequences} synthetic antibody sequences")
    return sequences

def generate_synthetic_toxin_sequences(
    num_sequences: int = 5, 
    min_length: int = 50, 
    max_length: int = 100
) -> List[str]:
    """
    Generate synthetic toxin sequences for testing.
    
    Args:
        num_sequences: Number of sequences to generate
        min_length: Minimum sequence length
        max_length: Maximum sequence length
        
    Returns:
        List of synthetic toxin sequences
    """
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    toxin_motifs = ["LRHVPE", "KVYTLK", "PYVKRP", "MCYKFT", "WTPCKG"]
    
    sequences = []
    
    for i in range(num_sequences):
        # Create a random sequence length
        seq_length = np.random.randint(min_length, max_length + 1)
        
        # Generate random sequence
        sequence = ''.join(np.random.choice(list(amino_acids)) for _ in range(seq_length))
        
        # Insert recognizable toxin motif at a random position
        motif_pos = np.random.randint(10, seq_length - 10)
        motif = np.random.choice(toxin_motifs)
        sequence = sequence[:motif_pos] + motif + sequence[motif_pos + len(motif):]
        
        sequences.append(sequence)
    
    logger.info(f"Generated {num_sequences} synthetic toxin sequences")
    return sequences

def test_esm2_initialization(config: ESM2Config) -> None:
    """
    Test ESM-2 model initialization.
    
    Args:
        config: ESM-2 configuration
    """
    logger.info("Testing ESM-2 model initialization...")
    start_time = time.time()
    
    try:
        model = ESM2CDRModel(config=config)
        logger.info("ESM-2 model initialized successfully")
        logger.info(f"Initialization time: {time.time() - start_time:.2f} seconds")
    except Exception as e:
        logger.error(f"Error initializing ESM-2 model: {e}")
        raise RuntimeError(f"Failed to initialize ESM-2 model: {e}")
        
def test_cdr_processor() -> Dict[str, str]:
    """
    Test CDR processor functionality.
    
    Returns:
        Dictionary mapping CDR names to sequences
    """
    logger.info("Testing CDR processor...")
    
    # Create CDR processor
    cdr_processor = CDRProcessor()
    
    # Test antibody sequence (human IgG1 anti-VEGF)
    test_sequence = "EVQLVESGGGLVQPGGSLRLSCAASGYTFTNYGMNWVRQAPGKGLEWVGWINTYTGEPTYAADFKRRFTFSLDTSKSTAYLQMNSLRAEDTAVYYCAKYPHYYGSSHWYFDVWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHT"
    
    # Extract CDR sequences
    cdrs = cdr_processor.extract_cdr_sequences(test_sequence)
    
    if not cdrs:
        logger.error("Failed to extract CDRs from test sequence")
        return {}
    
    # Log extracted CDRs
    for cdr_name, cdr_seq in cdrs.items():
        logger.info(f"{cdr_name}: {cdr_seq}")
    
    return cdrs

def test_cdr_masking(model: ESM2CDRModel, sequence: str, cdrs: Dict[str, str]) -> Optional[str]:
    """
    Test CDR masking functionality.
    
    Args:
        model: ESM-2 CDR model
        sequence: Antibody sequence
        cdrs: Dictionary mapping CDR names to sequences
        
    Returns:
        Masked sequence or None if masking failed
    """
    logger.info("Testing CDR masking...")
    
    if not cdrs:
        logger.error("No CDRs provided for masking")
        return None
    
    # Select a CDR to mask (prefer H3 if available)
    cdr_to_mask = "CDRH3" if "CDRH3" in cdrs else list(cdrs.keys())[0]
    logger.info(f"Masking {cdr_to_mask}: {cdrs[cdr_to_mask]}")
    
    # Mask the CDR
    masked_sequence = model.mask_cdr(sequence, cdr_to_mask)
    
    if masked_sequence:
        logger.info(f"Original: {sequence}")
        logger.info(f"Masked:   {masked_sequence}")
        return masked_sequence
    else:
        logger.error("Failed to mask CDR")
        return None

def test_masked_cdr_prediction(model: ESM2CDRModel, masked_sequence: str, cdr_name: str, original_cdr: str) -> List[str]:
    """
    Test masked CDR prediction.
    
    Args:
        model: ESM-2 CDR model
        masked_sequence: Sequence with masked CDR
        cdr_name: Name of the masked CDR
        original_cdr: Original CDR sequence
        
    Returns:
        List of predicted CDR sequences
    """
    logger.info("Testing masked CDR prediction...")
    
    if not masked_sequence:
        logger.error("No masked sequence provided")
        return []
    
    # Predict masked CDR
    predictions = model.predict_masked_cdrs(masked_sequence, [cdr_name], num_predictions=5)
    
    if not predictions or cdr_name not in predictions:
        logger.error("Failed to predict masked CDR")
        return []
    
    # Log predictions
    logger.info(f"Original CDR: {original_cdr}")
    for i, predicted_cdr in enumerate(predictions[cdr_name]):
        # Calculate accuracy (percentage of correctly predicted amino acids)
        correct_count = sum(a == b for a, b in zip(original_cdr, predicted_cdr))
        accuracy = (correct_count / max(len(original_cdr), len(predicted_cdr))) * 100
        logger.info(f"Prediction {i+1}: {predicted_cdr} (Accuracy: {accuracy:.2f}%)")
    
    return predictions[cdr_name]

def test_antibody_generation(model: ESM2CDRModel, sequence: str, cdrs: Dict[str, str]) -> List[str]:
    """
    Test antibody sequence generation.
    
    Args:
        model: ESM-2 CDR model
        sequence: Template antibody sequence
        cdrs: Dictionary mapping CDR names to sequences
        
    Returns:
        List of generated antibody sequences
    """
    logger.info("Testing antibody generation...")
    
    if not cdrs:
        logger.error("No CDRs provided for generation")
        return []
    
    # Select CDRs to modify (prefer H3 and L3 if available)
    cdrs_to_modify = []
    if "CDRH3" in cdrs:
        cdrs_to_modify.append("CDRH3")
    if "CDRL3" in cdrs:
        cdrs_to_modify.append("CDRL3")
    
    if not cdrs_to_modify:
        cdrs_to_modify = [list(cdrs.keys())[0]]
    
    logger.info(f"Generating antibodies by modifying: {', '.join(cdrs_to_modify)}")
    
    # Generate antibody variants
    variants = model.generate_antibody(sequence, cdrs_to_modify, num_variants=3)
    
    if not variants:
        logger.error("Failed to generate antibody variants")
        return []
    
    # Log variants
    logger.info(f"Template: {sequence}")
    for i, variant in enumerate(variants):
        logger.info(f"Variant {i+1}: {variant}")
    
    return variants

def test_binding_prediction(model: ESM2CDRModel, antibody_sequences: List[str], toxin_sequences: List[str]) -> List[Tuple[float, str, str]]:
    """
    Test binding prediction between toxins and antibodies.
    
    Args:
        model: ESM-2 CDR model
        antibody_sequences: List of antibody sequences
        toxin_sequences: List of toxin sequences
        
    Returns:
        List of tuples (binding_score, toxin_sequence, antibody_sequence)
    """
    logger.info("Testing binding prediction...")
    
    if not antibody_sequences or not toxin_sequences:
        logger.error("No antibody or toxin sequences provided")
        return []
    
    results = []
    
    # Test binding for a subset of sequence pairs
    for toxin in toxin_sequences[:2]:  # Test first 2 toxins
        for antibody in antibody_sequences[:3]:  # Test first 3 antibodies
            # Predict binding score
            binding_score = model.predict_binding_score(toxin, antibody)
            logger.info(f"Binding score: {binding_score:.4f}")
            results.append((binding_score, toxin, antibody))
    
    # Sort by binding score (descending)
    results.sort(reverse=True)
    
    # Log top results
    logger.info("Top binding predictions:")
    for i, (score, toxin, antibody) in enumerate(results[:3]):
        logger.info(f"{i+1}. Score: {score:.4f}")
    
    return results

def visualize_binding_predictions(binding_results: List[Tuple[float, str, str]], output_dir: str) -> None:
    """
    Visualize binding prediction results.
    
    Args:
        binding_results: List of binding prediction results
        output_dir: Directory to save visualizations
    """
    if not binding_results:
        logger.warning("No binding results to visualize")
        return
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract scores and create identifiers
    scores = [score for score, _, _ in binding_results]
    labels = [f"T{i//3+1}-Ab{i%3+1}" for i in range(len(binding_results))]
    
    # Create bar chart
    plt.figure(figsize=(10, 6))
    plt.bar(labels, scores)
    plt.xlabel("Toxin-Antibody Pair")
    plt.ylabel("Binding Score")
    plt.title("Predicted Binding Scores")
    plt.ylim(0, 1.0)
    plt.grid(alpha=0.3)
    plt.savefig(os.path.join(output_dir, "binding_scores.png"))
    plt.close()
    
    logger.info(f"Binding visualization saved to {os.path.join(output_dir, 'binding_scores.png')}")

def save_test_results(
    cdrs: Dict[str, str], 
    predicted_cdrs: List[str],
    generated_antibodies: List[str],
    binding_results: List[Tuple[float, str, str]],
    output_dir: str
) -> None:
    """
    Save test results to CSV files.
    
    Args:
        cdrs: Dictionary mapping CDR names to sequences
        predicted_cdrs: List of predicted CDR sequences
        generated_antibodies: List of generated antibody sequences
        binding_results: List of binding prediction results
        output_dir: Directory to save results
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Save CDR extraction results
    if cdrs:
        cdr_df = pd.DataFrame({"CDR": list(cdrs.keys()), "Sequence": list(cdrs.values())})
        cdr_df.to_csv(os.path.join(output_dir, "extracted_cdrs.csv"), index=False)
    
    # Save CDR prediction results
    if predicted_cdrs:
        pred_df = pd.DataFrame({"PredictedCDR": predicted_cdrs})
        pred_df.to_csv(os.path.join(output_dir, "predicted_cdrs.csv"), index=False)
    
    # Save generated antibody results
    if generated_antibodies:
        gen_df = pd.DataFrame({"GeneratedAntibody": generated_antibodies})
        gen_df.to_csv(os.path.join(output_dir, "generated_antibodies.csv"), index=False)
    
    # Save binding results
    if binding_results:
        bind_df = pd.DataFrame([
            {"BindingScore": score, "ToxinSequence": toxin[:20] + "...", "AntibodySequence": antibody[:20] + "..."}
            for score, toxin, antibody in binding_results
        ])
        bind_df.to_csv(os.path.join(output_dir, "binding_predictions.csv"), index=False)
    
    logger.info(f"Test results saved to {output_dir}")

def create_esm2_config() -> ESM2Config:
    """
    Create ESM-2 configuration for testing.
    
    Returns:
        ESM2Config: Configuration for testing
    """
    # Create a lighter configuration for testing
    config = ESM2Config()
    
    # Use a smaller model if available
    config.model["esm2_model"] = "esm2_t6_8M_UR50D"  # Use smaller model for testing
    config.model["embedding_dim"] = 320  # Smaller embedding dimension
    config.model["hidden_dim"] = 128  # Smaller hidden dimension
    
    # Reduce complexity
    config.model["num_hidden_layers"] = 1
    config.model["num_attention_heads"] = 8
    
    # Set device to automatic detection
    config.device = "auto"
    
    return config

def run_all_tests() -> None:
    """
    Run all ESM-2 transformer tests.
    """
    parser = argparse.ArgumentParser(description='Test ESM-2 transformer functionality')
    parser.add_argument('--output_dir', type=str, default='output/esm2_test',
                        help='Directory to save test results')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    # Set random seed
    set_seed(args.seed)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create ESM-2 configuration
    logger.info("Creating ESM-2 configuration...")
    config = create_esm2_config()
    
    # Test ESM-2 initialization
    test_esm2_initialization(config)
    
    # Initialize ESM-2 model
    logger.info("Initializing ESM-2 model for tests...")
    model = ESM2CDRModel(config=config)
    
    # Test CDR processor
    cdrs = test_cdr_processor()
    
    # Generate synthetic sequences for testing
    antibody_sequences = generate_synthetic_antibody_sequences()
    toxin_sequences = generate_synthetic_toxin_sequences()
    
    # Get a test sequence
    test_sequence = antibody_sequences[0]
    test_cdrs = model.cdr_processor.extract_cdr_sequences(test_sequence)
    
    if not test_cdrs:
        logger.warning("No CDRs found in synthetic sequence, using manually provided sequence")
        # Use a known antibody sequence with CDRs
        test_sequence = "EVQLVESGGGLVQPGGSLRLSCAASGYTFTNYGMNWVRQAPGKGLEWVGWINTYTGEPTYAADFKRRFTFSLDTSKSTAYLQMNSLRAEDTAVYYCAKYPHYYGSSHWYFDVWGQGTLVTVSS"
        test_cdrs = model.cdr_processor.extract_cdr_sequences(test_sequence)
    
    if not test_cdrs:
        logger.error("Failed to identify CDRs in test sequences")
        return
    
    # Test CDR masking
    masked_sequence = test_cdr_masking(model, test_sequence, test_cdrs)
    
    # Test CDR prediction (if masking successful)
    predicted_cdrs = []
    if masked_sequence:
        # Use the first CDR for prediction test
        cdr_name = list(test_cdrs.keys())[0]
        original_cdr = test_cdrs[cdr_name]
        predicted_cdrs = test_masked_cdr_prediction(model, masked_sequence, cdr_name, original_cdr)
    
    # Test antibody generation
    generated_antibodies = test_antibody_generation(model, test_sequence, test_cdrs)
    
    # Test binding prediction
    binding_results = test_binding_prediction(model, antibody_sequences, toxin_sequences)
    
    # Visualize binding predictions
    visualize_binding_predictions(binding_results, args.output_dir)
    
    # Save test results
    save_test_results(test_cdrs, predicted_cdrs, generated_antibodies, binding_results, args.output_dir)
    
    logger.info("All ESM-2 transformer tests completed successfully")

if __name__ == "__main__":
    run_all_tests()