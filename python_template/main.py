#!/usr/bin/env python
# main.py

import os
import sys
import logging
import argparse
import yaml
from typing import Dict, List, Any, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger("Phytovenomics.Main")

def setup_parser() -> argparse.ArgumentParser:
    """Set up command line argument parser."""
    parser = argparse.ArgumentParser(
        description="Phytovenomics ML Platform with ESM-2 Transformer"
    )
    
    # Main operation modes
    parser.add_argument(
        "--mode", 
        choices=["predict", "generate", "train", "demo"], 
        default="demo",
        help="Operation mode: predict CDRs, generate new CDRs, train model, or run demo"
    )
    
    # Model configuration
    parser.add_argument(
        "--config", 
        default="config/esm2_config.yaml", 
        help="Path to ESM-2 configuration file"
    )
    parser.add_argument(
        "--checkpoint", 
        help="Path to model checkpoint for inference or continued training"
    )
    
    # Input options
    parser.add_argument(
        "--sequence", 
        help="Input antibody sequence (for predict/generate modes)"
    )
    parser.add_argument(
        "--input_file", 
        help="File containing antibody sequences (one per line)"
    )
    parser.add_argument(
        "--cdr", 
        help="CDR to predict/generate (e.g., 'CDR-H3', 'CDR-L1')"
    )
    
    # Output options
    parser.add_argument(
        "--output_dir", 
        default="output", 
        help="Directory for output files"
    )
    parser.add_argument(
        "--output_prefix", 
        default="esm2_output", 
        help="Prefix for output files"
    )
    
    # Training options
    parser.add_argument(
        "--train_file", 
        help="Training data file (CSV/FASTA)"
    )
    parser.add_argument(
        "--val_file", 
        help="Validation data file (CSV/FASTA)"
    )
    parser.add_argument(
        "--epochs", 
        type=int, 
        default=10, 
        help="Number of training epochs"
    )
    
    # Generation options
    parser.add_argument(
        "--num_sequences", 
        type=int, 
        default=10, 
        help="Number of sequences to generate"
    )
    parser.add_argument(
        "--temperature", 
        type=float, 
        help="Temperature for sampling (overrides config)"
    )
    
    # Misc options
    parser.add_argument(
        "--seed", 
        type=int, 
        default=42, 
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--verbose", 
        action="store_true", 
        help="Enable verbose output"
    )
    
    return parser

def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from YAML file.
    
    Args:
        config_path: Path to the configuration file
        
    Returns:
        Configuration dictionary
    """
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        logger.info(f"Loaded configuration from {config_path}")
        return config
    except Exception as e:
        logger.error(f"Failed to load configuration: {e}")
        sys.exit(1)

def load_sequences(input_source: str) -> List[str]:
    """
    Load antibody sequences from a file or a single sequence.
    
    Args:
        input_source: Path to sequence file or a single sequence
        
    Returns:
        List of sequences
    """
    sequences = []
    
    if os.path.isfile(input_source):
        # Determine file type and load accordingly
        if input_source.endswith(('.fasta', '.fa')):
            from Bio import SeqIO
            for record in SeqIO.parse(input_source, "fasta"):
                sequences.append(str(record.seq))
        elif input_source.endswith('.csv'):
            import pandas as pd
            df = pd.read_csv(input_source)
            if 'sequence' in df.columns:
                sequences = df['sequence'].tolist()
            else:
                # Try to find a column that might contain sequences
                for col in df.columns:
                    if df[col].apply(lambda x: isinstance(x, str) and len(x) > 20).all():
                        sequences = df[col].tolist()
                        logger.info(f"Using column '{col}' as sequence data")
                        break
        else:
            # Assume text file with one sequence per line
            with open(input_source, 'r') as f:
                sequences = [line.strip() for line in f if line.strip()]
    else:
        # Assume input is a single sequence
        sequences = [input_source]
    
    return sequences

def predict_cdrs(model, sequences: List[str], cdr_names: Optional[List[str]] = None) -> Dict[str, List[Dict[str, Any]]]:
    """
    Predict masked CDRs in antibody sequences.
    
    Args:
        model: ESM2CDRModel instance
        sequences: List of antibody sequences
        cdr_names: Optional list of CDR names to predict
        
    Returns:
        Dictionary of prediction results
    """
    results = {}
    
    for i, sequence in enumerate(sequences):
        logger.info(f"Processing sequence {i+1}/{len(sequences)}")
        
        # Identify chain type and CDRs
        chain_type = model.cdr_processor.identify_chain_type(sequence)
        cdrs = model.cdr_processor.identify_cdrs(sequence)
        
        if not cdrs:
            logger.warning(f"No CDRs identified in sequence {i+1}")
            continue
        
        # If no specific CDRs are provided, predict all
        cdrs_to_predict = cdr_names if cdr_names else [name for _, _, name in cdrs]
        
        seq_results = []
        for cdr_name in cdrs_to_predict:
            # Check if this CDR exists in the sequence
            cdr_info = next((cdr for cdr in cdrs if cdr[2] == cdr_name), None)
            if not cdr_info:
                logger.warning(f"{cdr_name} not found in sequence {i+1}")
                continue
                
            start, end, _ = cdr_info
            cdr_sequence = sequence[start:end+1]
            
            # Mask and predict the CDR
            masked_seq, masked_positions = model.cdr_processor.mask_specific_cdrs(sequence, [cdr_name])
            predictions = model.predict_masked_residues(masked_seq, masked_positions)
            
            # Calculate accuracy
            correct = sum(1 for i, (pos, pred) in enumerate(zip(masked_positions, predictions)) 
                         if pred == sequence[pos])
            accuracy = 100 * correct / len(masked_positions) if masked_positions else 0
            
            # Store results
            seq_results.append({
                'cdr_name': cdr_name,
                'original_cdr': cdr_sequence,
                'predicted_cdr': ''.join(predictions),
                'positions': (start, end),
                'accuracy': accuracy,
                'length': len(cdr_sequence)
            })
        
        results[i] = seq_results
    
    return results

def generate_cdrs(model, sequences: List[str], num_variants: int = 1, 
                temperature: Optional[float] = None) -> Dict[str, List[Dict[str, Any]]]:
    """
    Generate new CDR sequences for antibodies.
    
    Args:
        model: ESM2CDRModel instance
        sequences: List of antibody sequences
        num_variants: Number of variants to generate per sequence
        temperature: Optional temperature override for sampling
        
    Returns:
        Dictionary of generation results
    """
    # Override temperature if provided
    if temperature is not None:
        original_temp = model.temperature
        model.temperature = temperature
    
    results = {}
    
    for i, sequence in enumerate(sequences):
        logger.info(f"Processing sequence {i+1}/{len(sequences)}")
        
        # Identify chain type
        chain_type = model.cdr_processor.identify_chain_type(sequence)
        
        # Set up CDR lengths based on chain type
        if chain_type == "heavy":
            cdr_names = ["CDR-H1", "CDR-H2", "CDR-H3"]
            base_lengths = {"CDR-H1": 8, "CDR-H2": 10, "CDR-H3": 15}
        elif chain_type == "light":
            cdr_names = ["CDR-L1", "CDR-L2", "CDR-L3"]
            base_lengths = {"CDR-L1": 11, "CDR-L2": 7, "CDR-L3": 9}
        else:
            logger.warning(f"Unknown chain type for sequence {i+1}")
            continue
        
        # Extract original CDRs
        original_cdrs = model.cdr_processor.extract_cdr_sequences(sequence)
        
        # Generate variants
        seq_results = []
        for variant in range(num_variants):
            # Vary CDR lengths slightly for each variant
            import random
            cdr_lengths = {}
            for cdr_name, base_length in base_lengths.items():
                # Use original length if available, otherwise use base length with small variations
                if cdr_name in original_cdrs:
                    orig_length = len(original_cdrs[cdr_name])
                    # Allow +/-2 variation from original length
                    var_length = orig_length + random.randint(-2, 2)
                    # Ensure reasonable bounds
                    min_len, max_len = model.cdr_processor.cdr_length_ranges[cdr_name]
                    cdr_lengths[cdr_name] = max(min_len, min(max_len, var_length))
                else:
                    cdr_lengths[cdr_name] = base_length
            
            # Generate new CDRs
            new_cdrs = model.generate_cdr_sequences(sequence, cdr_lengths)
            
            # Calculate properties for original and new CDRs
            cdr_properties = {}
            for cdr_name in cdr_names:
                if cdr_name in original_cdrs and cdr_name in new_cdrs:
                    orig_cdr = original_cdrs[cdr_name]
                    new_cdr = new_cdrs[cdr_name]
                    
                    # Calculate hydrophobicity
                    orig_hydro = model.cdr_processor.calculate_cdr_hydrophobicity(orig_cdr)
                    new_hydro = model.cdr_processor.calculate_cdr_hydrophobicity(new_cdr)
                    
                    cdr_properties[cdr_name] = {
                        'original': orig_cdr,
                        'generated': new_cdr,
                        'original_length': len(orig_cdr),
                        'generated_length': len(new_cdr),
                        'original_hydrophobicity': orig_hydro,
                        'generated_hydrophobicity': new_hydro,
                    }
            
            # Graft CDRs to get full sequence
            new_sequence = model.cdr_processor.graft_cdrs(sequence, new_cdrs, chain_type)
            
            # Store results
            seq_results.append({
                'variant': variant + 1,
                'original_sequence': sequence,
                'generated_sequence': new_sequence,
                'cdrs': cdr_properties
            })
        
        results[i] = seq_results
    
    # Restore original temperature if changed
    if temperature is not None:
        model.temperature = original_temp
    
    return results

def save_results(results: Dict, output_dir: str, prefix: str, mode: str):
    """
    Save results to output files.
    
    Args:
        results: Results dictionary
        output_dir: Output directory
        prefix: Output filename prefix
        mode: Operation mode
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Save detailed results as JSON
    import json
    json_path = os.path.join(output_dir, f"{prefix}_{mode}_detailed.json")
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    logger.info(f"Detailed results saved to {json_path}")
    
    # Save summary as CSV
    import pandas as pd
    
    if mode == "predict":
        # Create a list to hold prediction data
        predictions = []
        for seq_idx, seq_results in results.items():
            for result in seq_results:
                predictions.append({
                    'sequence_id': seq_idx,
                    'cdr_name': result['cdr_name'],
                    'original_cdr': result['original_cdr'],
                    'predicted_cdr': result['predicted_cdr'],
                    'accuracy': result['accuracy'],
                    'length': result['length']
                })
        
        if predictions:
            df = pd.DataFrame(predictions)
            csv_path = os.path.join(output_dir, f"{prefix}_predictions.csv")
            df.to_csv(csv_path, index=False)
            logger.info(f"Prediction summary saved to {csv_path}")
    
    elif mode == "generate":
        # Create a list to hold generation data
        generations = []
        for seq_idx, seq_results in results.items():
            for result in seq_results:
                for cdr_name, props in result['cdrs'].items():
                    generations.append({
                        'sequence_id': seq_idx,
                        'variant': result['variant'],
                        'cdr_name': cdr_name,
                        'original_cdr': props['original'],
                        'generated_cdr': props['generated'],
                        'original_length': props['original_length'],
                        'generated_length': props['generated_length'],
                        'original_hydrophobicity': props['original_hydrophobicity'],
                        'generated_hydrophobicity': props['generated_hydrophobicity']
                    })
        
        if generations:
            df = pd.DataFrame(generations)
            csv_path = os.path.join(output_dir, f"{prefix}_generations.csv")
            df.to_csv(csv_path, index=False)
            logger.info(f"Generation summary saved to {csv_path}")
            
            # Also save a FASTA file with generated sequences
            fasta_path = os.path.join(output_dir, f"{prefix}_generated_sequences.fasta")
            with open(fasta_path, 'w') as f:
                for seq_idx, seq_results in results.items():
                    for result in seq_results:
                        header = f">generated_seq_{seq_idx}_variant_{result['variant']}"
                        f.write(f"{header}\n{result['generated_sequence']}\n")
            logger.info(f"Generated sequences saved to {fasta_path}")

def train_model(model, train_file: str, val_file: Optional[str] = None, epochs: int = 10):
    """
    Train the ESM-2 model on antibody data.
    
    Args:
        model: ESM2CDRModel instance
        train_file: Path to training data file
        val_file: Optional path to validation data file
        epochs: Number of training epochs
    """
    from torch.utils.data import DataLoader, random_split
    from antibody_design.esm2_cdr_masking import AntibodyESMDataset
    
    # Load sequences
    train_sequences = load_sequences(train_file)
    logger.info(f"Loaded {len(train_sequences)} training sequences")
    
    val_sequences = []
    if val_file:
        val_sequences = load_sequences(val_file)
        logger.info(f"Loaded {len(val_sequences)} validation sequences")
    
    # Create datasets
    train_dataset = AntibodyESMDataset(
        train_sequences,
        model.cdr_processor,
        model.tokenizer,
        max_length=model.max_length,
        cdr_mask_prob=model.cdr_mask_prob,
        framework_mask_prob=model.framework_mask_prob
    )
    
    # If no validation file provided, split training data
    if not val_file:
        train_size = int(0.9 * len(train_dataset))
        val_size = len(train_dataset) - train_size
        train_dataset, val_dataset = random_split(train_dataset, [train_size, val_size])
    else:
        val_dataset = AntibodyESMDataset(
            val_sequences,
            model.cdr_processor,
            model.tokenizer,
            max_length=model.max_length,
            cdr_mask_prob=model.cdr_mask_prob,
            framework_mask_prob=model.framework_mask_prob
        )
    
    # Create data loaders
    train_dataloader = DataLoader(
        train_dataset,
        batch_size=model.batch_size,
        shuffle=True
    )
    
    val_dataloader = DataLoader(
        val_dataset,
        batch_size=model.batch_size
    )
    
    logger.info(f"Training on {len(train_dataset)} samples, validating on {len(val_dataset)} samples")
    logger.info(f"Training for {epochs} epochs")
    
    # Train the model
    model.train(train_dataloader, val_dataloader, epochs=epochs)
    
    logger.info("Training completed")

def run_demo(model, args):
    """
    Run demonstration of the ESM-2 model capabilities.
    
    Args:
        model: ESM2CDRModel instance
        args: Command line arguments
    """
    from antibody_design.esm2_cdr_demo import (
        analyze_framework_cdrs,
        predict_masked_cdrs,
        generate_new_cdrs
    )
    
    # If no sequence provided, use example sequences
    if args.sequence:
        sequence = args.sequence
    else:
        # Example heavy chain sequence
        heavy_chain = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRGDYYRYGMDVWGQGTTVTVSS"
        # Example light chain sequence
        light_chain = "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIK"
        
        sequence = heavy_chain
        logger.info("Using example heavy chain sequence")
    
    # Run demo functions
    analyze_framework_cdrs(model, sequence)
    
    if model.cdr_processor.identify_chain_type(sequence) == "heavy":
        predict_masked_cdrs(model, sequence, ["CDR-H3"])
    else:
        predict_masked_cdrs(model, sequence, ["CDR-L3"])
    
    generate_new_cdrs(model, sequence)

def main():
    """Main function to run the ESM-2 transformer functionality"""
    # Parse command line arguments
    parser = setup_parser()
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Load configuration
    config = load_config(args.config)
    
    # Import ESM-2 model
    try:
        from antibody_design.esm2_cdr_masking import load_model_from_config
        from utils.training_utils import set_seed
    except ImportError:
        logger.error("Failed to import required modules. Make sure you have installed all dependencies.")
        sys.exit(1)
    
    # Set random seed
    set_seed(args.seed)
    logger.info(f"Random seed set to {args.seed}")
    
    # Override config with command line arguments
    if args.temperature:
        config['generation']['temperature'] = args.temperature
    
    # Load model
    try:
        model = load_model_from_config(args.config)
        logger.info("ESM-2 model loaded successfully")
    except Exception as e:
        logger.error(f"Failed to load ESM-2 model: {e}")
        sys.exit(1)
    
    # Load checkpoint if provided
    if args.checkpoint:
        try:
            model.load_checkpoint(args.checkpoint)
            logger.info(f"Loaded checkpoint from {args.checkpoint}")
        except Exception as e:
            logger.error(f"Failed to load checkpoint: {e}")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run the specified mode
    if args.mode == "predict":
        # Check if sequence or input file is provided
        if not args.sequence and not args.input_file:
            logger.error("Either --sequence or --input_file must be provided for predict mode")
            sys.exit(1)
            
        # Load sequences
        sequences = []
        if args.sequence:
            sequences = [args.sequence]
        elif args.input_file:
            sequences = load_sequences(args.input_file)
            
        logger.info(f"Predicting CDRs for {len(sequences)} sequences")
        
        # Determine which CDRs to predict
        cdr_names = None
        if args.cdr:
            cdr_names = [args.cdr]
            
        # Perform prediction
        results = predict_cdrs(model, sequences, cdr_names)
        
        # Save results
        save_results(results, args.output_dir, args.output_prefix, "predict")
        
    elif args.mode == "generate":
        # Check if sequence or input file is provided
        if not args.sequence and not args.input_file:
            logger.error("Either --sequence or --input_file must be provided for generate mode")
            sys.exit(1)
            
        # Load sequences
        sequences = []
        if args.sequence:
            sequences = [args.sequence]
        elif args.input_file:
            sequences = load_sequences(args.input_file)
            
        logger.info(f"Generating CDRs for {len(sequences)} sequences")
        
        # Perform generation
        results = generate_cdrs(model, sequences, num_variants=args.num_sequences, 
                               temperature=args.temperature)
        
        # Save results
        save_results(results, args.output_dir, args.output_prefix, "generate")
        
    elif args.mode == "train":
        # Check if training file is provided
        if not args.train_file:
            logger.error("--train_file must be provided for train mode")
            sys.exit(1)
            
        # Train model
        train_model(model, args.train_file, args.val_file, args.epochs)
        
    elif args.mode == "demo":
        # Run demonstration
        run_demo(model, args)
    
    logger.info(f"ESM-2 {args.mode} completed successfully")

if __name__ == "__main__":
    main()