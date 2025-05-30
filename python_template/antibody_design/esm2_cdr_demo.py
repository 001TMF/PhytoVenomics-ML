#!/usr/bin/env python
# antibody_design/esm2_cdr_demo.py

import os
import sys
import logging
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import torch
from pathlib import Path

# Add parent directory to path if running as script
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

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
logger = logging.getLogger("Phytovenomics.ESM2Demo")

def create_output_directory(base_dir="output/esm2_demo"):
    """
    Create output directory for demo results.
    
    Args:
        base_dir: Base directory for outputs
        
    Returns:
        Path to the created directory
    """
    timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(base_dir, f"demo_{timestamp}")
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Created output directory: {output_dir}")
    return output_dir

def load_sample_data():
    """
    Load sample data for demonstration.
    
    Returns:
        Dictionary containing sample antibody and toxin sequences
    """
    # Example antibody sequences (human IgG1)
    antibodies = {
        "anti_vegf": "EVQLVESGGGLVQPGGSLRLSCAASGYTFTNYGMNWVRQAPGKGLEWVGWINTYTGEPTYAADFKRRFTFSLDTSKSTAYLQMNSLRAEDTAVYYCAKYPHYYGSSHWYFDVWGQGTLVTVSS",
        "anti_tnf": "EVQLVESGGGLVQPGGSLRLSCAASGYDFTHYGMNWVRQAPGKGLEWVGWINTYTGEPTYAADFKRRFTFSLDTSKSTAYLQMNSLRAEDTAVYYCAKYPYYYGTSHWYFDVWGQGTLVTVSS",
        "anti_il6": "EVQLVESGGGLVQPGGSLRLSCAASGFDFSRYWMSWVRQAPGKGLEWIGEINPDSSTINYAPSLKDKFIISRDNAKNSLYLQMNSLRAEDTAVYYCARWGGDGFYAMDYWGQGTLVTVSS"
    }
    
    # Example toxin sequences (snake venom toxins)
    toxins = {
        "cobra_neurotoxin": "LECHNQQSSQTPTTTGCSGGETNCYKKRWRDHRGYRTERGCGCPSVKNGIEINCCTTDRCNN",
        "viper_hemotoxin": "VVGGDECNINEHRFLVALYDGVTGVAYIGTSMCNRNTLDRPYNLGLYGEDGPDVLSYKSRRVNPL",
        "mamba_dendrotoxin": "RPRWIGGEGGPGWRQMERYVELPIVLNKRKMCNRDEPCSSIQIRWVGRCHICVFACPGPKRKIQVK"
    }
    
    return {"antibodies": antibodies, "toxins": toxins}

def demonstrate_cdr_identification(model):
    """
    Demonstrate CDR identification in antibody sequences.
    
    Args:
        model: ESM2CDRModel instance
        
    Returns:
        Dictionary of identified CDRs per antibody
    """
    logger.info("Demonstrating CDR identification...")
    
    # Load sample data
    data = load_sample_data()
    antibodies = data["antibodies"]
    
    # Extract CDRs from each antibody sequence
    results = {}
    for name, sequence in antibodies.items():
        logger.info(f"Processing {name} antibody sequence")
        cdrs = model.cdr_processor.extract_cdr_sequences(sequence)
        
        if cdrs:
            logger.info(f"Identified {len(cdrs)} CDR regions:")
            for cdr_name, cdr_seq in cdrs.items():
                logger.info(f"  {cdr_name}: {cdr_seq}")
            results[name] = cdrs
        else:
            logger.warning(f"No CDRs identified in {name}")
    
    return results

def demonstrate_cdr_masking(model, cdrs_by_antibody):
    """
    Demonstrate CDR masking in antibody sequences.
    
    Args:
        model: ESM2CDRModel instance
        cdrs_by_antibody: Dictionary of identified CDRs per antibody
        
    Returns:
        Dictionary of masked sequences
    """
    logger.info("Demonstrating CDR masking...")
    
    # Load sample data
    data = load_sample_data()
    antibodies = data["antibodies"]
    
    masked_sequences = {}
    
    for name, sequence in antibodies.items():
        if name not in cdrs_by_antibody or not cdrs_by_antibody[name]:
            logger.warning(f"No CDRs available for {name}, skipping")
            continue
            
        cdrs = cdrs_by_antibody[name]
        
        # Mask all CDRs
        logger.info(f"Masking all CDRs in {name}")
        all_masked = model.mask_cdrs(sequence, list(cdrs.keys()))
        if all_masked:
            masked_sequences[f"{name}_all_masked"] = all_masked
            logger.info(f"All CDRs masked: {all_masked}")
        
        # Mask individual CDRs
        for cdr_name in cdrs:
            logger.info(f"Masking {cdr_name} in {name}")
            masked = model.mask_cdr(sequence, cdr_name)
            if masked:
                masked_sequences[f"{name}_{cdr_name}_masked"] = masked
                logger.info(f"{cdr_name} masked: {masked}")
    
    return masked_sequences

def demonstrate_cdr_prediction(model, masked_sequences):
    """
    Demonstrate CDR prediction from masked sequences.
    
    Args:
        model: ESM2CDRModel instance
        masked_sequences: Dictionary of masked antibody sequences
        
    Returns:
        Dictionary of predicted CDRs
    """
    logger.info("Demonstrating CDR prediction...")
    
    predictions = {}
    
    for name, sequence in masked_sequences.items():
        # Skip sequences with all CDRs masked
        if "_all_masked" in name:
            continue
            
        logger.info(f"Predicting CDR for {name}")
        
        # Extract the CDR name from the masked sequence name
        parts = name.split("_")
        if len(parts) >= 3:
            cdr_name = parts[-2]
            if not cdr_name.startswith("CDR"):
                cdr_name = "CDR" + cdr_name
                
            # Predict masked CDR
            cdr_predictions = model.predict_masked_cdrs(sequence, [cdr_name], num_predictions=3)
            
            if cdr_predictions and cdr_name in cdr_predictions:
                logger.info(f"Predicted sequences for {cdr_name}:")
                for i, pred in enumerate(cdr_predictions[cdr_name]):
                    logger.info(f"  Prediction {i+1}: {pred}")
                predictions[name] = cdr_predictions[cdr_name]
            else:
                logger.warning(f"Failed to predict {cdr_name} for {name}")
    
    return predictions

def demonstrate_antibody_generation(model, num_variants=3):
    """
    Demonstrate antibody generation by modifying CDRs.
    
    Args:
        model: ESM2CDRModel instance
        num_variants: Number of variants to generate
        
    Returns:
        Dictionary of generated antibody variants
    """
    logger.info("Demonstrating antibody generation...")
    
    # Load sample data
    data = load_sample_data()
    antibodies = data["antibodies"]
    
    generated_antibodies = {}
    
    for name, sequence in antibodies.items():
        logger.info(f"Generating variants for {name}")
        
        # Extract CDRs
        cdrs = model.cdr_processor.extract_cdr_sequences(sequence)
        
        if not cdrs:
            logger.warning(f"No CDRs identified in {name}, skipping")
            continue
            
        # Define CDRs to modify (prefer H3 and L3 if available)
        cdrs_to_modify = []
        if "CDRH3" in cdrs:
            cdrs_to_modify.append("CDRH3")
        if "CDRL3" in cdrs:
            cdrs_to_modify.append("CDRL3")
        
        if not cdrs_to_modify:
            # If H3/L3 not found, use first available CDR
            cdrs_to_modify = [list(cdrs.keys())[0]]
            
        logger.info(f"Modifying CDRs: {', '.join(cdrs_to_modify)}")
        
        # Generate variants
        variants = model.generate_antibody(sequence, cdrs_to_modify, num_variants=num_variants)
        
        if variants:
            logger.info(f"Generated {len(variants)} variants")
            for i, variant in enumerate(variants):
                logger.info(f"Variant {i+1}: {variant[:50]}...")  # Show beginning of sequence
            generated_antibodies[name] = variants
        else:
            logger.warning(f"Failed to generate variants for {name}")
    
    return generated_antibodies

def demonstrate_binding_prediction(model):
    """
    Demonstrate binding prediction between toxins and antibodies.
    
    Args:
        model: ESM2CDRModel instance
        
    Returns:
        DataFrame of binding predictions
    """
    logger.info("Demonstrating binding prediction...")
    
    # Load sample data
    data = load_sample_data()
    antibodies = data["antibodies"]
    toxins = data["toxins"]
    
    # Create dataframe to store results
    results = []
    
    for toxin_name, toxin_seq in toxins.items():
        for ab_name, ab_seq in antibodies.items():
            logger.info(f"Predicting binding between {toxin_name} and {ab_name}")
            
            # Predict binding score
            binding_score = model.predict_binding_score(toxin_seq, ab_seq)
            
            results.append({
                "Toxin": toxin_name,
                "Antibody": ab_name,
                "BindingScore": binding_score
            })
            
            logger.info(f"Binding score: {binding_score:.4f}")
    
    # Convert to DataFrame
    binding_df = pd.DataFrame(results)
    
    return binding_df

def visualize_binding_scores(binding_df, output_dir):
    """
    Visualize binding scores between toxins and antibodies.
    
    Args:
        binding_df: DataFrame of binding predictions
        output_dir: Directory to save visualizations
    """
    logger.info("Visualizing binding scores...")
    
    # Create pivot table for heatmap
    pivot_df = binding_df.pivot(index="Toxin", columns="Antibody", values="BindingScore")
    
    # Create heatmap
    plt.figure(figsize=(10, 8))
    plt.imshow(pivot_df.values, cmap="YlOrRd")
    plt.colorbar(label="Binding Score")
    plt.xticks(range(len(pivot_df.columns)), pivot_df.columns, rotation=45, ha="right")
    plt.yticks(range(len(pivot_df.index)), pivot_df.index)
    plt.title("Predicted Binding Scores between Toxins and Antibodies")
    
    # Save figure
    heatmap_path = os.path.join(output_dir, "binding_heatmap.png")
    plt.tight_layout()
    plt.savefig(heatmap_path)
    plt.close()
    
    logger.info(f"Saved binding heatmap to {heatmap_path}")
    
    # Create bar chart for each toxin
    for toxin in binding_df["Toxin"].unique():
        toxin_data = binding_df[binding_df["Toxin"] == toxin]
        
        plt.figure(figsize=(10, 6))
        bars = plt.bar(toxin_data["Antibody"], toxin_data["BindingScore"])
        plt.ylim(0, 1.0)
        plt.xlabel("Antibody")
        plt.ylabel("Binding Score")
        plt.title(f"Binding Scores for {toxin}")
        plt.grid(alpha=0.3)
        
        # Add value labels on top of bars
        for bar in bars:
            height = bar.get_height()
            plt.annotate(f'{height:.2f}',
                         xy=(bar.get_x() + bar.get_width() / 2, height),
                         xytext=(0, 3),
                         textcoords="offset points",
                         ha='center', va='bottom')
        
        # Save figure
        chart_path = os.path.join(output_dir, f"binding_scores_{toxin}.png")
        plt.tight_layout()
        plt.savefig(chart_path)
        plt.close()
        
        logger.info(f"Saved binding chart for {toxin} to {chart_path}")

def save_results_to_csv(output_dir, cdrs, masked_sequences, predictions, generated_antibodies, binding_df):
    """
    Save demonstration results to CSV files.
    
    Args:
        output_dir: Directory to save results
        cdrs: Dictionary of identified CDRs
        masked_sequences: Dictionary of masked sequences
        predictions: Dictionary of predicted CDRs
        generated_antibodies: Dictionary of generated antibody variants
        binding_df: DataFrame of binding predictions
    """
    logger.info("Saving results to CSV files...")
    
    # Save CDRs
    cdr_rows = []
    for ab_name, cdr_dict in cdrs.items():
        for cdr_name, cdr_seq in cdr_dict.items():
            cdr_rows.append({
                "Antibody": ab_name,
                "CDR": cdr_name,
                "Sequence": cdr_seq
            })
    
    if cdr_rows:
        cdr_df = pd.DataFrame(cdr_rows)
        cdr_df.to_csv(os.path.join(output_dir, "identified_cdrs.csv"), index=False)
    
    # Save masked sequences
    masked_df = pd.DataFrame([
        {"Name": name, "MaskedSequence": seq} 
        for name, seq in masked_sequences.items()
    ])
    if not masked_df.empty:
        masked_df.to_csv(os.path.join(output_dir, "masked_sequences.csv"), index=False)
    
    # Save predictions
    pred_rows = []
    for name, preds in predictions.items():
        for i, pred in enumerate(preds):
            pred_rows.append({
                "Name": name,
                "PredictionRank": i+1,
                "PredictedSequence": pred
            })
    
    if pred_rows:
        pred_df = pd.DataFrame(pred_rows)
        pred_df.to_csv(os.path.join(output_dir, "cdr_predictions.csv"), index=False)
    
    # Save generated antibodies
    gen_rows = []
    for ab_name, variants in generated_antibodies.items():
        for i, variant in enumerate(variants):
            gen_rows.append({
                "TemplateAntibody": ab_name,
                "VariantNumber": i+1,
                "Sequence": variant
            })
    
    if gen_rows:
        gen_df = pd.DataFrame(gen_rows)
        gen_df.to_csv(os.path.join(output_dir, "generated_antibodies.csv"), index=False)
    
    # Save binding predictions
    if not binding_df.empty:
        binding_df.to_csv(os.path.join(output_dir, "binding_predictions.csv"), index=False)
    
    logger.info(f"Results saved to {output_dir}")

def create_demo_summary(output_dir, cdrs, masked_sequences, predictions, generated_antibodies, binding_df):
    """
    Create a summary report of the demonstration.
    
    Args:
        output_dir: Directory to save summary
        cdrs: Dictionary of identified CDRs
        masked_sequences: Dictionary of masked sequences
        predictions: Dictionary of predicted CDRs
        generated_antibodies: Dictionary of generated antibody variants
        binding_df: DataFrame of binding predictions
    """
    logger.info("Creating demo summary report...")
    
    summary_path = os.path.join(output_dir, "demo_summary.md")
    
    with open(summary_path, 'w') as f:
        f.write("# ESM-2 Transformer Demonstration Summary\n\n")
        f.write(f"Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # CDR Identification
        f.write("## 1. CDR Identification\n\n")
        f.write(f"Identified CDRs in {len(cdrs)} antibody sequences.\n\n")
        for ab_name, cdr_dict in cdrs.items():
            f.write(f"### {ab_name}\n\n")
            f.write("| CDR | Sequence |\n")
            f.write("|-----|----------|\n")
            for cdr_name, cdr_seq in cdr_dict.items():
                f.write(f"| {cdr_name} | {cdr_seq} |\n")
            f.write("\n")
        
        # CDR Masking
        f.write("## 2. CDR Masking\n\n")
        f.write(f"Created {len(masked_sequences)} masked sequences.\n\n")
        f.write("Examples:\n\n")
        for name, seq in list(masked_sequences.items())[:3]:
            f.write(f"- **{name}**: {seq}\n")
        f.write("\n")
        
        # CDR Prediction
        f.write("## 3. CDR Prediction\n\n")
        f.write(f"Generated predictions for {len(predictions)} masked CDRs.\n\n")
        for name, preds in predictions.items():
            f.write(f"### {name}\n\n")
            f.write("| Rank | Predicted Sequence |\n")
            f.write("|------|-------------------|\n")
            for i, pred in enumerate(preds):
                f.write(f"| {i+1} | {pred} |\n")
            f.write("\n")
        
        # Antibody Generation
        f.write("## 4. Antibody Generation\n\n")
        f.write(f"Generated variants for {len(generated_antibodies)} antibody templates.\n\n")
        for ab_name, variants in generated_antibodies.items():
            f.write(f"### {ab_name}\n\n")
            f.write("| Variant | Sequence |\n")
            f.write("|---------|----------|\n")
            for i, variant in enumerate(variants):
                # Show first 50 characters followed by ellipsis
                f.write(f"| {i+1} | {variant[:50]}... |\n")
            f.write("\n")
        
        # Binding Prediction
        f.write("## 5. Binding Prediction\n\n")
        f.write(f"Predicted binding scores between {binding_df['Toxin'].nunique()} toxins and {binding_df['Antibody'].nunique()} antibodies.\n\n")
        f.write("### Top Binding Pairs\n\n")
        
        # Get top 5 binding pairs
        top_binding = binding_df.nlargest(5, "BindingScore")
        f.write("| Toxin | Antibody | Binding Score |\n")
        f.write("|-------|----------|---------------|\n")
        for _, row in top_binding.iterrows():
            f.write(f"| {row['Toxin']} | {row['Antibody']} | {row['BindingScore']:.4f} |\n")
        f.write("\n")
        
        # References to files
        f.write("## Output Files\n\n")
        f.write("The following output files were generated:\n\n")
        f.write("- `identified_cdrs.csv`: List of all identified CDRs\n")
        f.write("- `masked_sequences.csv`: Sequences with masked CDRs\n")
        f.write("- `cdr_predictions.csv`: Predictions for masked CDRs\n")
        f.write("- `generated_antibodies.csv`: Generated antibody variants\n")
        f.write("- `binding_predictions.csv`: Binding scores between toxins and antibodies\n")
        f.write("- Visualizations:\n")
        f.write("  - `binding_heatmap.png`: Heatmap of binding scores\n")
        f.write("  - `binding_scores_*.png`: Bar charts of binding scores per toxin\n")
    
    logger.info(f"Demo summary report saved to {summary_path}")
    return summary_path

def run_demonstration(config_path=None):
    """
    Run the full ESM-2 transformer demonstration.
    
    Args:
        config_path: Path to configuration file (optional)
    """
    # Create output directory
    output_dir = create_output_directory()
    
    # Create or load configuration
    if config_path and os.path.exists(config_path):
        logger.info(f"Loading configuration from {config_path}")
        config = ESM2Config.from_yaml(config_path)
    else:
        logger.info("Using default configuration")
        config = ESM2Config()
        
        # Use smaller model for demonstration to reduce memory requirements
        config.model["esm2_model"] = "esm2_t6_8M_UR50D"  # Use smaller model
        config.model["embedding_dim"] = 320  # Smaller embedding dimension
    
    # Set random seed
    set_seed(config.training["seed"])
    
    try:
        # Initialize model
        logger.info("Initializing ESM-2 model...")
        model = ESM2CDRModel(config=config)
        
        # Run demonstration components
        cdrs = demonstrate_cdr_identification(model)
        masked_sequences = demonstrate_cdr_masking(model, cdrs)
        predictions = demonstrate_cdr_prediction(model, masked_sequences)
        generated_antibodies = demonstrate_antibody_generation(model)
        binding_df = demonstrate_binding_prediction(model)
        
        # Visualize results
        visualize_binding_scores(binding_df, output_dir)
        
        # Save results
        save_results_to_csv(output_dir, cdrs, masked_sequences, predictions, generated_antibodies, binding_df)
        
        # Create summary
        summary_path = create_demo_summary(output_dir, cdrs, masked_sequences, predictions, generated_antibodies, binding_df)
        
        logger.info("Demonstration completed successfully!")
        logger.info(f"Results saved to {output_dir}")
        logger.info(f"Summary report: {summary_path}")
        
        return output_dir, summary_path
        
    except Exception as e:
        logger.error(f"Demonstration failed: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Demonstrate ESM-2 CDR masking and antibody design')
    parser.add_argument('--config', type=str, help='Path to configuration file')
    
    args = parser.parse_args()
    run_demonstration(config_path=args.config)