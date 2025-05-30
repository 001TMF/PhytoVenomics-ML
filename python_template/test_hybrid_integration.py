#!/usr/bin/env python
# test_hybrid_integration.py

import os
import sys
import logging
import argparse
import yaml
import torch
import tempfile
import json
import random
import time
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("Phytovenomics.HybridIntegrationTest")

# Add the parent directory to the Python path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import Phytovenomics modules
from antibody_design.evolutionary_search import EvolutionarySearch, Individual, RosettaFoldClient
from antibody_design.evolutionary_hybrid_integration import HybridEvolutionarySearch, StructureBasedEvaluator

# Import structure prediction modules
from src.structure_prediction.hybrid_prediction import HybridPredictionPipeline

# Import utility modules
from utils.validation_utils import ValidationMetrics
from utils.visualization import plot_antibody_structures
from utils.training_utils import set_seed

class HybridIntegrationTest:
    """
    Test suite to validate the integration between the evolutionary search module
    and the hybrid prediction pipeline in the Phytovenomics ML platform.
    """
    
    def __init__(self, config_path: Optional[str] = None, output_dir: Optional[str] = None, use_mock: bool = True):
        """
        Initialize the hybrid integration test.
        
        Args:
            config_path: Path to configuration file
            output_dir: Directory for test outputs
            use_mock: Whether to use mock implementations for heavy computations
        """
        # Set random seed for reproducibility
        set_seed(42)
        
        # Set up output directory
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = Path("test_output/hybrid_integration")
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Set mock mode
        self.use_mock = use_mock
        os.environ["USE_MOCK_PREDICTION"] = "True" if use_mock else "False"
        
        # Load configuration
        self.config = self._load_config(config_path)
        
        # Initialize test data
        self.test_epitope = None
        self.test_antibodies = None
        
        # Metrics for test results
        self.metrics = ValidationMetrics()
        
        logger.info(f"Hybrid Integration Test initialized with output to {self.output_dir}")
        
    def _load_config(self, config_path: Optional[str] = None) -> Dict[str, Any]:
        """
        Load configuration from file or use default test configuration.
        
        Args:
            config_path: Path to configuration file
            
        Returns:
            Configuration dictionary
        """
        if config_path and os.path.exists(config_path):
            try:
                with open(config_path, "r") as f:
                    config = yaml.safe_load(f)
                logger.info(f"Loaded configuration from {config_path}")
                return config
            except Exception as e:
                logger.warning(f"Failed to load configuration from {config_path}: {e}")
        
        # Use default test configuration
        logger.info("Using default test configuration")
        return {
            "output_dir": str(self.output_dir),
            "use_mock": self.use_mock,
            "evolutionary_search": {
                "population_size": 8,
                "generations": 3,
                "mutation_rate": 0.1,
                "crossover_rate": 0.7,
                "elite_size": 2,
            },
            "hybrid_prediction": {
                "enabled": True,
                "methods": ["igfold", "esmfold"],
                "confidence_threshold": 0.7,
                "batch_size": 2
            },
            "rosetta": {
                "rosetta_path": "/opt/rosetta/",  # Default path, will be ignored in mock mode
                "scripts_path": "/opt/rosetta/scripts/",
                "use_gpu": False,  # Use CPU for testing
                "use_mock": self.use_mock
            }
        }
    
    def create_test_data(self) -> None:
        """
        Create synthetic test data for the integration test.
        """
        # Create test epitope
        self.test_epitope = {
            "id": "TEST_EPITOPE_001",
            "sequence": "RVQPTESIVRFPNITNLCPFGEVFNATR",
            "start": 10,
            "end": 37,
            "score": 0.92,
            "toxin_id": "TOX001",
            "toxin_family": "three_finger_toxin",
            "source_species": "Naja naja",
            "prediction_method": "ML-enhanced"
        }
        
        # Create test antibodies (initial population)
        self.test_antibodies = []
        
        # Heavy and light chain framework templates
        heavy_framework = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS---WVRQAPGKGLEWVSAINSGS---YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR------WGQGTLVTVSS"
        light_framework = "DIQMTQSPSSLSASVGDRVTITCRASQ------WYQQKPGKAPKLLIY-------GVPSRFSGSGSGTDFTLTISSLQPEDFATYYC--------FGQGTKVEIK"
        
        # CDR sequences (variable regions)
        cdr_h1_options = ["SYGMH", "SYNMH", "SYGMT", "SYWMS"]
        cdr_h2_options = ["GISGSGGSTY", "GIIPIFGTA", "GIIPMFGTA", "GIDPSFGTA"]
        cdr_h3_options = ["DYGDYGDYFDYWG", "DGGDYAMDYWG", "REGGVFDYWG", "AKGWFDYWG", "DYGDYFDYWG"]
        
        cdr_l1_options = ["SISSYLN", "SITNYLA", "SITSYLN", "SIQSYLN"]
        cdr_l2_options = ["AASSLQS", "AASRLQS", "AASNLQS", "AASRLHS"]
        cdr_l3_options = ["QQSYSTPPT", "QQSYSTSYT", "QQSYSIPPT", "QQSYSIPYT"]
        
        # Generate 8 diverse antibodies
        for i in range(8):
            # Select CDR sequences
            cdr_h1 = random.choice(cdr_h1_options)
            cdr_h2 = random.choice(cdr_h2_options)
            cdr_h3 = random.choice(cdr_h3_options)
            cdr_l1 = random.choice(cdr_l1_options)
            cdr_l2 = random.choice(cdr_l2_options)
            cdr_l3 = random.choice(cdr_l3_options)
            
            # Insert CDRs into framework
            heavy_chain = heavy_framework.replace("---", cdr_h1, 1)
            heavy_chain = heavy_chain.replace("---", cdr_h2, 1)
            heavy_chain = heavy_chain.replace("------", cdr_h3, 1)
            
            light_chain = light_framework.replace("------", cdr_l1, 1)
            light_chain = light_framework.replace("-------", cdr_l2, 1)
            light_chain = light_framework.replace("--------", cdr_l3, 1)
            
            # Create individual
            ind = Individual(
                heavy_chain=heavy_chain,
                light_chain=light_chain,
                cdr_sequences={
                    "H1": cdr_h1,
                    "H2": cdr_h2,
                    "H3": cdr_h3,
                    "L1": cdr_l1,
                    "L2": cdr_l2,
                    "L3": cdr_l3
                },
                name=f"TEST_AB_{i+1:03d}"
            )
            self.test_antibodies.append(ind)
            
        logger.info(f"Created test data: {len(self.test_antibodies)} antibodies and 1 epitope")
        
        # Save test data to file
        test_data = {
            "epitope": self.test_epitope,
            "antibodies": [
                {
                    "name": ab.name,
                    "heavy_chain": ab.heavy_chain,
                    "light_chain": ab.light_chain,
                    "cdr_sequences": ab.cdr_sequences
                }
                for ab in self.test_antibodies
            ]
        }
        
        with open(self.output_dir / "test_data.json", "w") as f:
            json.dump(test_data, f, indent=2)
    
    def test_structure_based_evaluator(self) -> bool:
        """
        Test the structure-based evaluator component.
        
        Returns:
            Whether the test passed
        """
        logger.info("Testing StructureBasedEvaluator...")
        
        try:
            # Create test directory
            evaluator_dir = self.output_dir / "structure_evaluator"
            os.makedirs(evaluator_dir, exist_ok=True)
            
            # Create structure-based evaluator
            evaluator = StructureBasedEvaluator(
                config=self.config,
                use_mock=self.use_mock,
                cache_dir=str(evaluator_dir / "cache")
            )
            
            # Test structure prediction
            test_antibody = self.test_antibodies[0]
            logger.info(f"Predicting structure for {test_antibody.name}...")
            structure_result = evaluator.predict_structure(
                individual=test_antibody,
                output_dir=str(evaluator_dir)
            )
            
            # Check result
            if not structure_result.get("success", False):
                logger.error(f"Structure prediction failed: {structure_result.get('error', 'Unknown error')}")
                self.metrics.add_metric("structure_prediction_success", 0)
                return False
                
            logger.info(f"Structure predicted successfully using {structure_result.get('method_used', 'unknown')}")
            self.metrics.add_metric("structure_prediction_success", 1)
            self.metrics.add_metric("prediction_time", structure_result.get("prediction_time", 0))
            
            # Test binding evaluation
            binding_energy = evaluator.evaluate_binding(
                individual=test_antibody,
                target_epitope=self.test_epitope,
                structure_result=structure_result
            )
            logger.info(f"Binding energy: {binding_energy:.2f}")
            self.metrics.add_metric("binding_energy", binding_energy)
            
            # Test stability evaluation
            stability = evaluator.evaluate_stability(
                individual=test_antibody,
                structure_result=structure_result
            )
            logger.info(f"Stability: {stability:.2f}")
            self.metrics.add_metric("stability", stability)
            
            # Test developability evaluation
            developability = evaluator.evaluate_developability(
                individual=test_antibody,
                structure_result=structure_result
            )
            logger.info(f"Developability: {developability:.2f}")
            self.metrics.add_metric("developability", developability)
            
            return True
            
        except Exception as e:
            logger.error(f"Error in structure-based evaluator test: {str(e)}")
            return False
    
    def test_hybrid_search_integration(self) -> bool:
        """
        Test the integration between evolutionary search and hybrid prediction pipeline.
        
        Returns:
            Whether the test passed
        """
        logger.info("Testing HybridEvolutionarySearch integration...")
        
        try:
            # Create test directory
            hybrid_dir = self.output_dir / "hybrid_search"
            os.makedirs(hybrid_dir, exist_ok=True)
            
            # Update config to point to the test directory
            self.config["output_dir"] = str(hybrid_dir)
            
            # Save config to file for use with HybridEvolutionarySearch
            config_path = hybrid_dir / "hybrid_config.yaml"
            with open(config_path, "w") as f:
                yaml.dump(self.config, f)
                
            # Create hybrid evolutionary search
            start_time = time.time()
            hybrid_search = HybridEvolutionarySearch(str(config_path))
            
            # Run optimization
            logger.info(f"Starting hybrid evolutionary search with {len(self.test_antibodies)} initial antibodies...")
            optimized_antibodies = hybrid_search.optimize_antibody(
                initial_population=self.test_antibodies,
                target_epitope=self.test_epitope
            )
            
            # Record timing
            end_time = time.time()
            optimization_time = end_time - start_time
            self.metrics.add_metric("optimization_time", optimization_time)
            
            # Check results
            if not optimized_antibodies:
                logger.error("No optimized antibodies returned")
                return False
                
            logger.info(f"Hybrid evolutionary search completed in {optimization_time:.2f}s")
            logger.info(f"Optimized {len(optimized_antibodies)} antibodies")
            
            # Analyze results
            best_ab = optimized_antibodies[0]  # Should be sorted by fitness
            logger.info(f"Best antibody: {best_ab.name}")
            logger.info(f"  - Binding energy: {best_ab.binding_energy:.2f}")
            logger.info(f"  - Stability: {best_ab.stability:.2f}")
            logger.info(f"  - Developability: {best_ab.developability:.2f}")
            
            # Save optimization results
            results = {
                "optimization_time": optimization_time,
                "population_size": len(self.test_antibodies),
                "optimized_size": len(optimized_antibodies),
                "best_antibody": {
                    "name": best_ab.name,
                    "binding_energy": best_ab.binding_energy,
                    "stability": best_ab.stability,
                    "developability": best_ab.developability,
                    "heavy_chain": best_ab.heavy_chain,
                    "light_chain": best_ab.light_chain,
                    "cdr_sequences": best_ab.cdr_sequences
                },
                "all_antibodies": [
                    {
                        "name": ab.name,
                        "binding_energy": ab.binding_energy,
                        "stability": ab.stability,
                        "developability": ab.developability
                    }
                    for ab in optimized_antibodies
                ]
            }
            
            with open(hybrid_dir / "optimization_results.json", "w") as f:
                json.dump(results, f, indent=2)
                
            # Record best antibody metrics
            self.metrics.add_metric("best_binding_energy", best_ab.binding_energy)
            self.metrics.add_metric("best_stability", best_ab.stability)
            self.metrics.add_metric("best_developability", best_ab.developability)
            
            return True
            
        except Exception as e:
            logger.error(f"Error in hybrid search integration test: {str(e)}")
            return False
    
    def test_hybrid_prediction_pipeline(self) -> bool:
        """
        Test the hybrid prediction pipeline component directly.
        
        Returns:
            Whether the test passed
        """
        logger.info("Testing HybridPredictionPipeline...")
        
        try:
            # Create test directory
            pipeline_dir = self.output_dir / "prediction_pipeline"
            os.makedirs(pipeline_dir, exist_ok=True)
            
            # Initialize hybrid prediction pipeline
            pipeline = HybridPredictionPipeline(
                use_cache=True,
                cache_dir=str(pipeline_dir / "cache")
            )
            
            # Test antibody sequence
            test_ab = self.test_antibodies[0]
            
            # Predict structure for heavy chain
            start_time = time.time()
            heavy_result = pipeline.predict_antibody_structure(
                sequence=test_ab.heavy_chain,
                chain_type="heavy"
            )
            heavy_time = time.time() - start_time
            
            # Predict structure for light chain
            start_time = time.time()
            light_result = pipeline.predict_antibody_structure(
                sequence=test_ab.light_chain,
                chain_type="light"
            )
            light_time = time.time() - start_time
            
            # Record metrics
            self.metrics.add_metric("heavy_chain_prediction_time", heavy_time)
            self.metrics.add_metric("light_chain_prediction_time", light_time)
            self.metrics.add_metric("heavy_chain_confidence", heavy_result.confidence)
            self.metrics.add_metric("light_chain_confidence", light_result.confidence)
            
            # Analyze structure quality
            heavy_quality = pipeline.analyze_structure_quality(heavy_result)
            light_quality = pipeline.analyze_structure_quality(light_result)
            
            logger.info(f"Heavy chain prediction: {heavy_time:.2f}s, confidence: {heavy_result.confidence:.2f}")
            logger.info(f"Light chain prediction: {light_time:.2f}s, confidence: {light_result.confidence:.2f}")
            
            # Visualize structures if possible
            try:
                heavy_viz = pipeline.visualize_structure(
                    heavy_result, 
                    output_path=str(pipeline_dir / "heavy_chain_structure.png"),
                    visualization_type="matplotlib"
                )
                
                light_viz = pipeline.visualize_structure(
                    light_result, 
                    output_path=str(pipeline_dir / "light_chain_structure.png"),
                    visualization_type="matplotlib"
                )
                
                logger.info("Structure visualizations created")
                
            except Exception as viz_e:
                logger.warning(f"Could not create structure visualizations: {viz_e}")
            
            # Save results
            results = {
                "heavy_chain": {
                    "prediction_time": heavy_time,
                    "confidence": heavy_result.confidence,
                    "model_name": heavy_result.model_name,
                    "quality_metrics": heavy_quality
                },
                "light_chain": {
                    "prediction_time": light_time,
                    "confidence": light_result.confidence,
                    "model_name": light_result.model_name,
                    "quality_metrics": light_quality
                }
            }
            
            with open(pipeline_dir / "prediction_results.json", "w") as f:
                json.dump(results, f, indent=2)
                
            return True
            
        except Exception as e:
            logger.error(f"Error in hybrid prediction pipeline test: {str(e)}")
            return False
            
    def save_metrics(self) -> None:
        """
        Save test metrics to file.
        """
        metrics_json = self.metrics.export_to_json()
        
        with open(self.output_dir / "test_metrics.json", "w") as f:
            f.write(metrics_json)
            
        # Create plots
        self.metrics.plot_metrics(
            ["binding_energy", "best_binding_energy"],
            output_path=str(self.output_dir / "binding_energy_metrics.png"),
            title="Binding Energy Metrics"
        )
        
        self.metrics.plot_metrics(
            ["stability", "best_stability"],
            output_path=str(self.output_dir / "stability_metrics.png"),
            title="Stability Metrics"
        )
        
        self.metrics.plot_metrics(
            ["heavy_chain_prediction_time", "light_chain_prediction_time", "optimization_time"],
            output_path=str(self.output_dir / "timing_metrics.png"),
            title="Timing Metrics"
        )
        
        logger.info(f"Test metrics saved to {self.output_dir / 'test_metrics.json'}")
    
    def run_all_tests(self) -> Dict[str, bool]:
        """
        Run all hybrid integration tests.
        
        Returns:
            Dictionary of test results
        """
        logger.info("Starting hybrid integration tests...")
        
        # Create test data
        self.create_test_data()
        
        # Run tests
        results = {}
        
        # Test the structure-based evaluator
        results["structure_evaluator"] = self.test_structure_based_evaluator()
        
        # Test the hybrid prediction pipeline
        results["prediction_pipeline"] = self.test_hybrid_prediction_pipeline()
        
        # Test the hybrid evolutionary search integration
        results["hybrid_search"] = self.test_hybrid_search_integration()
        
        # Save metrics
        self.save_metrics()
        
        # Print summary
        logger.info("\n----- Hybrid Integration Test Summary -----")
        for test, passed in results.items():
            status = "PASSED" if passed else "FAILED"
            logger.info(f"{test}: {status}")
        logger.info("----------------------------------------\n")
        
        # Create summary file
        with open(self.output_dir / "test_summary.txt", "w") as f:
            f.write("Hybrid Integration Test Summary\n")
            f.write("==============================\n\n")
            for test, passed in results.items():
                status = "PASSED" if passed else "FAILED"
                f.write(f"{test}: {status}\n")
            f.write("\n")
            
        return results

def main():
    """Main function to run hybrid integration tests."""
    parser = argparse.ArgumentParser(description="Phytovenomics Hybrid Integration Test")
    
    parser.add_argument("--config", help="Path to configuration file")
    parser.add_argument("--output-dir", default="test_output/hybrid_integration", help="Output directory")
    parser.add_argument("--no-mock", action="store_true", help="Disable mock mode (use real computations)")
    parser.add_argument("--test", choices=["all", "structure", "pipeline", "hybrid"], default="all", help="Specific test to run")
    parser.add_argument("--verbose", "-v", action="count", default=0, help="Increase verbosity")
    
    args = parser.parse_args()
    
    # Set verbosity
    if args.verbose >= 2:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.verbose == 1:
        logging.getLogger().setLevel(logging.INFO)
    
    # Create test instance
    tester = HybridIntegrationTest(
        config_path=args.config,
        output_dir=args.output_dir,
        use_mock=not args.no_mock
    )
    
    # Create test data
    tester.create_test_data()
    
    # Run specified test or all tests
    if args.test == "all":
        results = tester.run_all_tests()
        success = all(results.values())
    elif args.test == "structure":
        success = tester.test_structure_based_evaluator()
    elif args.test == "pipeline":
        success = tester.test_hybrid_prediction_pipeline()
    elif args.test == "hybrid":
        success = tester.test_hybrid_search_integration()
    
    # Save metrics
    tester.save_metrics()
    
    # Return exit code
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())