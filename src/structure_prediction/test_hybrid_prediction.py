#!/usr/bin/env python
# test_hybrid_prediction.py

import os
import sys
import time
import logging
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from typing import Dict, List, Optional, Tuple
import tempfile
import json
import argparse
from concurrent.futures import ThreadPoolExecutor

# Add the project root to Python path to import modules correctly
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

# Import the hybrid prediction pipeline
from src.structure_prediction.hybrid_prediction import (
    HybridPredictionPipeline,
    PredictionResult,
    ESMFoldInterface,
    IgFoldInterface
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("Phytovenomics.TestHybridPrediction")

# Constants
DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data"
ML_READY_DATA_DIR = DATA_DIR / "ml_ready"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
ANTIBODY_DATA_FILE = ML_READY_DATA_DIR / "antibodies_final_train.csv"

def load_antibody_sequences(limit: int = 10) -> pd.DataFrame:
    """
    Load antibody sequences from the training dataset
    
    Args:
        limit: Maximum number of sequences to load
    
    Returns:
        DataFrame containing antibody data
    """
    logger.info(f"Loading antibody data from {ANTIBODY_DATA_FILE}")
    try:
        df = pd.read_csv(ANTIBODY_DATA_FILE)
        logger.info(f"Loaded {len(df)} antibody sequences")
        
        # Filter for heavy chains for this test
        heavy_chains = df[df['chain_type'] == 'Heavy Chain'].head(limit)
        logger.info(f"Selected {len(heavy_chains)} heavy chain sequences for testing")
        
        return heavy_chains
    except Exception as e:
        logger.error(f"Error loading antibody data: {e}")
        raise

def test_model_availability():
    """Test if ESMFold and IgFold models are available"""
    logger.info("Testing model availability...")
    
    # Test ESMFold
    try:
        esm_interface = ESMFoldInterface()
        if esm_interface.is_available():
            logger.info("✅ ESMFold model is available")
        else:
            logger.warning("⚠️ ESMFold model is not available")
    except Exception as e:
        logger.error(f"❌ Error initializing ESMFold: {e}")
    
    # Test IgFold
    try:
        igfold_interface = IgFoldInterface()
        if igfold_interface.is_available():
            logger.info("✅ IgFold model is available")
        else:
            logger.warning("⚠️ IgFold model is not available")
    except Exception as e:
        logger.error(f"❌ Error initializing IgFold: {e}")

def run_single_prediction_test(sequence: str, description: str):
    """
    Run prediction test for a single antibody sequence
    
    Args:
        sequence: Antibody amino acid sequence
        description: Description of the antibody
    """
    logger.info(f"Testing prediction for: {description[:50]}...")
    
    pipeline = HybridPredictionPipeline(use_cache=True)
    
    # Test with ESMFold
    try:
        logger.info("Predicting structure with ESMFold...")
        start_time = time.time()
        esm_result = pipeline.predict_antibody_structure(
            sequence=sequence, 
            prefer_model="esmfold"
        )
        esm_time = time.time() - start_time
        logger.info(f"ESMFold prediction completed in {esm_time:.2f} seconds")
        logger.info(f"ESMFold confidence: {esm_result.confidence:.4f}")
    except Exception as e:
        logger.error(f"ESMFold prediction failed: {e}")
        esm_result = None
    
    # Test with IgFold
    try:
        logger.info("Predicting structure with IgFold...")
        start_time = time.time()
        igfold_result = pipeline.predict_antibody_structure(
            sequence=sequence, 
            prefer_model="igfold", 
            chain_type="heavy"
        )
        igfold_time = time.time() - start_time
        logger.info(f"IgFold prediction completed in {igfold_time:.2f} seconds")
        logger.info(f"IgFold confidence: {igfold_result.confidence:.4f}")
    except Exception as e:
        logger.error(f"IgFold prediction failed: {e}")
        igfold_result = None
    
    # Test ensemble approach
    try:
        logger.info("Predicting structure with hybrid ensemble...")
        start_time = time.time()
        ensemble_result = pipeline.predict_antibody_structure(
            sequence=sequence, 
            use_ensemble=True
        )
        ensemble_time = time.time() - start_time
        logger.info(f"Ensemble prediction completed in {ensemble_time:.2f} seconds")
        logger.info(f"Ensemble model used: {ensemble_result.model_name}")
        logger.info(f"Ensemble confidence: {ensemble_result.confidence:.4f}")
    except Exception as e:
        logger.error(f"Ensemble prediction failed: {e}")
        ensemble_result = None
    
    # Create results directory if it doesn't exist
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    # Save visualizations
    if esm_result and igfold_result:
        try:
            # Create unique identifier for this test
            test_id = description[:20].replace(" ", "_").replace(",", "")
            
            # Save ESMFold visualization
            esm_viz_path = str(RESULTS_DIR / f"{test_id}_esmfold.png")
            pipeline.visualize_structure(
                esm_result, 
                output_path=esm_viz_path,
                visualization_type="matplotlib"
            )
            
            # Save IgFold visualization
            igfold_viz_path = str(RESULTS_DIR / f"{test_id}_igfold.png")
            pipeline.visualize_structure(
                igfold_result, 
                output_path=igfold_viz_path,
                visualization_type="matplotlib", 
                highlight_cdrs=True
            )
            
            logger.info(f"Saved visualizations to {RESULTS_DIR}")
            
            # Return comparative results
            return {
                "description": description,
                "sequence_length": len(sequence),
                "esm_time": esm_time if esm_result else None,
                "igfold_time": igfold_time if igfold_result else None,
                "ensemble_time": ensemble_time if ensemble_result else None,
                "esm_confidence": esm_result.confidence if esm_result else None,
                "igfold_confidence": igfold_result.confidence if igfold_result else None,
                "ensemble_confidence": ensemble_result.confidence if ensemble_result else None,
                "ensemble_model": ensemble_result.model_name if ensemble_result else None
            }
        except Exception as e:
            logger.error(f"Error saving visualizations: {e}")
    
    return None

def benchmark_prediction_speed(sequences: List[str], chain_types: List[str], num_runs: int = 3):
    """
    Benchmark prediction speed for ESMFold and IgFold
    
    Args:
        sequences: List of antibody sequences to test
        chain_types: List of chain types for each sequence
        num_runs: Number of runs for each test to average results
    """
    logger.info(f"Benchmarking prediction speed with {len(sequences)} sequences...")
    
    pipeline = HybridPredictionPipeline(use_cache=False)  # Disable cache for benchmarking
    
    results = {
        "esm_times": [],
        "igfold_times": [],
        "sequence_lengths": []
    }
    
    for i, (sequence, chain_type) in enumerate(zip(sequences, chain_types)):
        logger.info(f"Testing sequence {i+1}/{len(sequences)} (length: {len(sequence)})")
        
        # Run ESMFold benchmark
        esm_times = []
        for run in range(num_runs):
            try:
                start_time = time.time()
                _ = pipeline.predict_antibody_structure(
                    sequence=sequence, 
                    prefer_model="esmfold"
                )
                esm_times.append(time.time() - start_time)
                logger.info(f"ESMFold run {run+1}/{num_runs}: {esm_times[-1]:.2f}s")
            except Exception as e:
                logger.error(f"ESMFold benchmark failed: {e}")
        
        # Run IgFold benchmark
        igfold_times = []
        for run in range(num_runs):
            try:
                start_time = time.time()
                _ = pipeline.predict_antibody_structure(
                    sequence=sequence, 
                    prefer_model="igfold",
                    chain_type=chain_type.lower()
                )
                igfold_times.append(time.time() - start_time)
                logger.info(f"IgFold run {run+1}/{num_runs}: {igfold_times[-1]:.2f}s")
            except Exception as e:
                logger.error(f"IgFold benchmark failed: {e}")
        
        # Record results
        if esm_times:
            results["esm_times"].append(sum(esm_times) / len(esm_times))
        if igfold_times:
            results["igfold_times"].append(sum(igfold_times) / len(igfold_times))
        results["sequence_lengths"].append(len(sequence))
    
    # Generate benchmark plots
    if results["esm_times"] and results["igfold_times"]:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        
        plt.figure(figsize=(10, 6))
        plt.scatter(results["sequence_lengths"], results["esm_times"], label="ESMFold", marker="o", alpha=0.7)
        plt.scatter(results["sequence_lengths"], results["igfold_times"], label="IgFold", marker="x", alpha=0.7)
        
        # Add trend lines
        if len(results["sequence_lengths"]) > 1:
            esm_z = np.polyfit(results["sequence_lengths"], results["esm_times"], 1)
            esm_p = np.poly1d(esm_z)
            plt.plot(results["sequence_lengths"], esm_p(results["sequence_lengths"]), "b--", alpha=0.5)
            
            igfold_z = np.polyfit(results["sequence_lengths"], results["igfold_times"], 1)
            igfold_p = np.poly1d(igfold_z)
            plt.plot(results["sequence_lengths"], igfold_p(results["sequence_lengths"]), "r--", alpha=0.5)
        
        plt.xlabel("Sequence Length")
        plt.ylabel("Prediction Time (s)")
        plt.title("Structure Prediction Time Comparison")
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        plt.savefig(str(RESULTS_DIR / "prediction_benchmark.png"), dpi=300, bbox_inches="tight")
        logger.info(f"Saved benchmark plot to {RESULTS_DIR / 'prediction_benchmark.png'}")
        
        # Calculate average speedup
        avg_speedup = sum(e/i for e, i in zip(results["esm_times"], results["igfold_times"])) / len(results["esm_times"])
        logger.info(f"Average speedup of IgFold over ESMFold: {avg_speedup:.2f}x")
        
        # Save benchmark results to JSON
        with open(RESULTS_DIR / "benchmark_results.json", "w") as f:
            json.dump({
                "sequence_lengths": results["sequence_lengths"],
                "esm_times": results["esm_times"],
                "igfold_times": results["igfold_times"],
                "average_speedup": avg_speedup
            }, f, indent=2)

def compare_cdr_accuracy():
    """
    Compare CDR region accuracy between ESMFold and IgFold
    This test requires ground truth structures which we don't have available,
    so this is a placeholder for future implementation.
    """
    logger.info("CDR accuracy comparison test is not implemented yet")
    logger.info("This would require ground truth structures with known CDR regions")
    logger.info("Future implementation would compare RMSD values for CDR regions")

def generate_test_report(test_results: List[Dict]):
    """
    Generate a test report from the results
    
    Args:
        test_results: List of test result dictionaries
    """
    if not test_results:
        logger.warning("No test results available for report generation")
        return
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    report_path = RESULTS_DIR / "test_report.md"
    
    with open(report_path, "w") as f:
        f.write("# Hybrid Prediction Pipeline Test Report\n\n")
        f.write(f"Test Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Model Availability\n\n")
        f.write("- ESMFold: Available\n")
        f.write("- IgFold: Available\n\n")
        
        f.write("## Performance Comparison\n\n")
        
        f.write("| Antibody | Length | ESMFold Time (s) | IgFold Time (s) | ESMFold Conf. | IgFold Conf. |\n")
        f.write("|----------|--------|------------------|-----------------|---------------|-------------|\n")
        
        for result in test_results:
            f.write(f"| {result['description'][:20]}... | {result['sequence_length']} | ")
            f.write(f"{result['esm_time']:.2f} | {result['igfold_time']:.2f} | ")
            f.write(f"{result['esm_confidence']:.4f} | {result['igfold_confidence']:.4f} |\n")
        
        # Calculate averages
        avg_esm_time = sum(r['esm_time'] for r in test_results if r['esm_time']) / len([r for r in test_results if r['esm_time']])
        avg_igfold_time = sum(r['igfold_time'] for r in test_results if r['igfold_time']) / len([r for r in test_results if r['igfold_time']])
        avg_esm_conf = sum(r['esm_confidence'] for r in test_results if r['esm_confidence']) / len([r for r in test_results if r['esm_confidence']])
        avg_igfold_conf = sum(r['igfold_confidence'] for r in test_results if r['igfold_confidence']) / len([r for r in test_results if r['igfold_confidence']])
        
        f.write("| **Average** | - | ")
        f.write(f"{avg_esm_time:.2f} | {avg_igfold_time:.2f} | ")
        f.write(f"{avg_esm_conf:.4f} | {avg_igfold_conf:.4f} |\n\n")
        
        f.write("## Speed Comparison\n\n")
        speedup = avg_esm_time / avg_igfold_time
        f.write(f"Average speedup of IgFold over ESMFold: **{speedup:.2f}x**\n\n")
        
        f.write("## Ensemble Model Selection\n\n")
        ensemble_choices = {}
        for result in test_results:
            model = result.get('ensemble_model', 'Unknown')
            ensemble_choices[model] = ensemble_choices.get(model, 0) + 1
            
        for model, count in ensemble_choices.items():
            f.write(f"- {model}: {count} selections ({count/len(test_results)*100:.1f}%)\n")
        
        f.write("\n## Conclusion\n\n")
        if speedup >= 10:
            f.write("IgFold demonstrates **significantly faster** structure prediction compared to ESMFold, ")
        elif speedup >= 2:
            f.write("IgFold demonstrates **moderately faster** structure prediction compared to ESMFold, ")
        else:
            f.write("IgFold demonstrates **comparable speed** to ESMFold, ")
            
        if avg_igfold_conf > avg_esm_conf:
            f.write("while also achieving higher average confidence scores. ")
        else:
            f.write("but with slightly lower average confidence scores. ")
            
        if "IgFold" in ensemble_choices and ensemble_choices.get("IgFold", 0) / len(test_results) > 0.7:
            f.write("The ensemble approach predominantly selects IgFold predictions, ")
            f.write("indicating its superior performance for antibody structure prediction.\n")
        else:
            f.write("The ensemble approach selects a mix of both models, ")
            f.write("indicating complementary strengths in different scenarios.\n")
    
    logger.info(f"Test report generated at {report_path}")
    return report_path

def main():
    """Main function for testing the hybrid prediction pipeline"""
    parser = argparse.ArgumentParser(description="Test the hybrid prediction pipeline")
    parser.add_argument("--num-sequences", type=int, default=3, help="Number of antibody sequences to test")
    parser.add_argument("--benchmark", action="store_true", help="Run performance benchmark")
    parser.add_argument("--quick", action="store_true", help="Run quick test with minimal sequences")
    args = parser.parse_args()
    
    # Create results directory
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    # Test model availability
    test_model_availability()
    
    # Load antibody sequences
    limit = 1 if args.quick else args.num_sequences
    antibodies = load_antibody_sequences(limit=limit)
    
    if antibodies.empty:
        logger.error("No antibody sequences available for testing")
        return
    
    # Run prediction tests
    test_results = []
    for _, row in antibodies.iterrows():
        result = run_single_prediction_test(
            sequence=row['sequence'],
            description=f"{row['description']} ({row['id']})"
        )
        if result:
            test_results.append(result)
    
    # Generate test report
    if test_results:
        report_path = generate_test_report(test_results)
        logger.info(f"Test report available at: {report_path}")
    else:
        logger.warning("No test results available, skipping report generation")
    
    # Run benchmark if requested
    if args.benchmark and not args.quick:
        benchmark_prediction_speed(
            sequences=antibodies['sequence'].tolist(),
            chain_types=["heavy"] * len(antibodies),
            num_runs=3
        )
    
    logger.info("Testing completed successfully")

if __name__ == "__main__":
    main()