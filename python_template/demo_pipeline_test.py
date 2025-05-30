#!/usr/bin/env python
# demo_pipeline_test.py

import os
import sys
import logging
import time
import json
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional
import matplotlib.pyplot as plt

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("Phytovenomics.DemoTest")

# Add src directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import components from antibody_design module
try:
    from antibody_design.evolutionary_search import EvolutionarySearch, Individual
    from antibody_design.evolutionary_hybrid_integration import HybridEvolutionarySearch
    from utils.validation_metrics import ValidationMetrics
    import utils.visualization as viz  # Import as module, not class
except ImportError as e:
    logger.error(f"Failed to import core modules: {e}")
    logger.error("Make sure you're running this script from the python_template directory")
    sys.exit(1)

# Import structure prediction module
try:
    from src.structure_prediction.hybrid_prediction import HybridPredictionPipeline
except ImportError as e:
    logger.error(f"Failed to import structure prediction module: {e}")
    logger.error("Creating mock implementation for testing")
    
    # Create mock class for testing
    class HybridPredictionPipeline:
        def __init__(self, use_cache=True, cache_dir=None):
            self.use_cache = use_cache
            self.cache_dir = cache_dir
            logger.warning("Using mock HybridPredictionPipeline")
            
        def predict_antibody_structure(self, sequence, chain_type="heavy", use_ensemble=True):
            from collections import namedtuple
            PredictionResult = namedtuple("PredictionResult", 
                                         ["pdb_data", "confidence", "prediction_time", 
                                          "model_name", "metadata", "sequence"])
            return PredictionResult(
                pdb_data="MOCK PDB DATA",
                confidence=0.85,
                prediction_time=1.5,
                model_name="mock_model",
                metadata={"mock": True},
                sequence=sequence
            )
            
        def analyze_structure_quality(self, result):
            return {"confidence": result.confidence, "quality": "medium"}
            
        def visualize_structure(self, result, output_path=None, visualization_type=None, highlight_cdrs=False):
            if output_path:
                with open(output_path, 'w') as f:
                    f.write("Mock structure visualization")
            return output_path or "mock_visualization"

# Import cocktail strategy module (with fallback mock)
try:
    from cocktail_strategy import CocktailOptimizer
except ImportError as e:
    logger.error(f"Failed to import cocktail strategy module: {e}")
    logger.error("Creating mock implementation for testing")
    
    # Create mock class for testing
    class CocktailOptimizer:
        def __init__(self, config=None):
            self.config = config or {}
            logger.warning("Using mock CocktailOptimizer")
            
        def design_optimal_cocktail(self, antibodies, toxins):
            covered_toxins = {toxin_id: 0.8 for toxin_id in toxins.keys()}
            return {
                "antibodies": [ab["id"] for ab in antibodies[:3]],
                "coverage": covered_toxins,
                "average_coverage": 0.8,
                "rationale": "Mock cocktail optimization"
            }


class PipelineDemoTest:
    """
    Demonstration test for the Phytovenomics pipeline integration.
    """
    
    def __init__(self, output_dir: str = "demo_test_results", use_mocks: bool = True):
        """
        Initialize the demo test.
        
        Args:
            output_dir: Directory for test outputs
            use_mocks: Whether to use mock implementations when components are unavailable
        """
        # Set up output directory
        self.output_dir = Path(output_dir)
        os.makedirs(self.output_dir, exist_ok=True)
        self.use_mocks = use_mocks
        
        # Configure environment for mocks if needed
        if use_mocks:
            os.environ["USE_MOCK_PREDICTION"] = "True"
            os.environ["USE_MOCK_EVALUATION"] = "True"
        
        # Initialize metrics tracker
        self.metrics = ValidationMetrics(str(self.output_dir))
        
        logger.info(f"Demo test initialized with output directory: {self.output_dir}")
        logger.info(f"Mock mode: {self.use_mocks}")
        
    def create_test_data(self) -> Dict[str, Any]:
        """
        Create synthetic test data for the pipeline demo.
        
        Returns:
            Dictionary containing test data
        """
        logger.info("Creating test data...")
        
        # Create test epitope
        epitope = {
            "id": "DEMO_EPITOPE_001",
            "sequence": "GFNCYFPLQSYGFQPTNGVGYQ",
            "start": 435,
            "end": 456,
            "score": 0.89,
            "toxin_id": "TOXIN001",
            "toxin_family": "three_finger_toxin",
            "source_species": "Naja naja"
        }
        
        # Create test toxins
        toxins = {
            "TOXIN001": {
                "id": "TOXIN001",
                "sequence": "LKCNKLVPLFYKTCPAGKNLCYKMFMVATPKVPVKRGCIDVCPKSSLLVKYVCCNTDRCN",
                "family": "three_finger_toxin",
                "source_species": "Naja naja",
                "function": "neurotoxic",
                "potency": 2.3
            },
            "TOXIN002": {
                "id": "TOXIN002",
                "sequence": "MKTLLLTLVVVTIVCLDLGYTRICYNHQSTTRATTKSCEENSCYKKYWRDHRGTIIERGCGCPTVKPGIKLSCCESEVCNN",
                "family": "phospholipase_a2",
                "source_species": "Bungarus caeruleus",
                "function": "neurotoxic",
                "potency": 1.8
            }
        }
        
        # Create initial antibody population
        # Define framework regions
        heavy_framework = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS---WVRQAPGKGLEWVSAINSGS---YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR------WGQGTLVTVSS"
        light_framework = "DIQMTQSPSSLSASVGDRVTITCRASQ------WYQQKPGKAPKLLIY-------GVPSRFSGSGSGTDFTLTISSLQPEDFATYYC--------FGQGTKVEIK"
        
        # Define CDR sequence options (simplified)
        cdr_h1_options = ["SYGMY", "SYWMS", "SYGMT", "SYNMH"]
        cdr_h2_options = ["TISSGGSYT", "GIDPSFGTA", "GIIPMFGTA", "GISGSGGSTY"]
        cdr_h3_options = ["AKGWFDYWG", "DYGDYFDYWG", "DGGDYAMDYWG", "REGGVFDYWG"]
        
        cdr_l1_options = ["SISSYLN", "SITNYLA", "SIQSYLN", "SITSYLN"]
        cdr_l2_options = ["AASSLQS", "AASRLQS", "AASNLQS", "AASRLHS"]
        cdr_l3_options = ["QQSYSTPPT", "QQSYSTSYT", "QQSYSIPPT", "QQSYSIPYT"]
        
        # Generate initial population
        initial_antibodies = []
        for i in range(6):
            # Select random CDR sequences for diversity
            import random
            cdr_h1 = random.choice(cdr_h1_options)
            cdr_h2 = random.choice(cdr_h2_options)
            cdr_h3 = random.choice(cdr_h3_options)
            cdr_l1 = random.choice(cdr_l1_options)
            cdr_l2 = random.choice(cdr_l2_options)
            cdr_l3 = random.choice(cdr_l3_options)
            
            # Insert CDRs into frameworks
            heavy_chain = heavy_framework.replace("---", cdr_h1, 1)
            heavy_chain = heavy_chain.replace("---", cdr_h2, 1)
            heavy_chain = heavy_chain.replace("------", cdr_h3, 1)
            
            light_chain = light_framework.replace("------", cdr_l1, 1)
            light_chain = light_framework.replace("-------", cdr_l2, 1)
            light_chain = light_framework.replace("--------", cdr_l3, 1)
            
            # Create Individual object
            initial_antibodies.append(Individual(
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
                name=f"DEMO_AB_{i+1:02d}"
            ))
        
        # Save test data to file
        test_data = {
            "epitope": epitope,
            "toxins": list(toxins.values()),
            "initial_antibodies": [
                {
                    "name": ab.name,
                    "heavy_chain": ab.heavy_chain,
                    "light_chain": ab.light_chain,
                    "cdr_sequences": ab.cdr_sequences
                }
                for ab in initial_antibodies
            ]
        }
        
        with open(self.output_dir / "test_data.json", "w") as f:
            json.dump(test_data, f, indent=2)
        
        logger.info(f"Created test data with {len(initial_antibodies)} antibodies and {len(toxins)} toxins")
        
        return {
            "epitope": epitope,
            "toxins": toxins,
            "initial_antibodies": initial_antibodies
        }
        
    def test_structure_prediction(self, antibody: Individual) -> Dict[str, Any]:
        """
        Test the structure prediction pipeline.
        
        Args:
            antibody: Antibody individual to predict structure for
            
        Returns:
            Dictionary with test results
        """
        logger.info(f"Testing structure prediction for {antibody.name}...")
        
        # Create structure prediction output directory
        structure_dir = self.output_dir / "structure_prediction"
        os.makedirs(structure_dir, exist_ok=True)
        
        try:
            # Initialize structure prediction pipeline
            structure_pipeline = HybridPredictionPipeline(
                use_cache=True,
                cache_dir=str(structure_dir / "cache")
            )
            
            # Time the prediction
            start_time = time.time()
            
            # Predict heavy chain structure
            heavy_result = structure_pipeline.predict_antibody_structure(
                sequence=antibody.heavy_chain,
                chain_type="heavy",
                use_ensemble=True
            )
            
            # Predict light chain structure
            light_result = structure_pipeline.predict_antibody_structure(
                sequence=antibody.light_chain,
                chain_type="light",
                use_ensemble=True
            )
            
            # Calculate total time
            total_time = time.time() - start_time
            
            # Analyze structure quality
            heavy_quality = structure_pipeline.analyze_structure_quality(heavy_result)
            light_quality = structure_pipeline.analyze_structure_quality(light_result)
            
            # Create visualizations if possible
            try:
                heavy_viz_path = str(structure_dir / f"{antibody.name}_heavy_chain.png")
                light_viz_path = str(structure_dir / f"{antibody.name}_light_chain.png")
                
                structure_pipeline.visualize_structure(
                    heavy_result,
                    output_path=heavy_viz_path,
                    visualization_type="matplotlib",
                    highlight_cdrs=True
                )
                
                structure_pipeline.visualize_structure(
                    light_result,
                    output_path=light_viz_path,
                    visualization_type="matplotlib",
                    highlight_cdrs=True
                )
                
                logger.info("Structure visualizations created successfully")
            except Exception as viz_err:
                logger.warning(f"Failed to create visualizations: {viz_err}")
            
            # Record metrics
            self.metrics.add_metric("heavy_chain_prediction_time", heavy_result.prediction_time)
            self.metrics.add_metric("light_chain_prediction_time", light_result.prediction_time)
            self.metrics.add_metric("total_prediction_time", total_time)
            self.metrics.add_metric("heavy_chain_confidence", heavy_result.confidence)
            self.metrics.add_metric("light_chain_confidence", light_result.confidence)
            
            # Create results dictionary
            results = {
                "success": True,
                "heavy_chain": {
                    "prediction_time": heavy_result.prediction_time,
                    "confidence": heavy_result.confidence,
                    "model_name": heavy_result.model_name,
                    "quality": heavy_quality
                },
                "light_chain": {
                    "prediction_time": light_result.prediction_time,
                    "confidence": light_result.confidence,
                    "model_name": light_result.model_name,
                    "quality": light_quality
                },
                "total_time": total_time
            }
            
            # Save results
            with open(structure_dir / f"{antibody.name}_prediction_results.json", "w") as f:
                json.dump(results, f, indent=2)
                
            logger.info(f"Structure prediction completed in {total_time:.2f}s")
            return results
            
        except Exception as e:
            logger.error(f"Structure prediction failed: {e}")
            return {"success": False, "error": str(e)}
    
    def test_hybrid_evolutionary_search(self, initial_antibodies: List[Individual], target_epitope: Dict) -> Dict[str, Any]:
        """
        Test the hybrid evolutionary search pipeline.
        
        Args:
            initial_antibodies: List of initial antibody designs
            target_epitope: Target epitope dictionary
            
        Returns:
            Dictionary with test results
        """
        logger.info("Testing hybrid evolutionary search...")
        
        # Create evolutionary search output directory
        evolution_dir = self.output_dir / "evolutionary_search"
        os.makedirs(evolution_dir, exist_ok=True)
        
        try:
            # Create config for evolutionary search
            config = {
                "population_size": len(initial_antibodies),
                "generations": 3,  # Keep small for demo
                "mutation_rate": 0.1,
                "crossover_rate": 0.7,
                "elite_size": 1,
                "output_dir": str(evolution_dir),
                "hybrid_prediction": {
                    "enabled": True,
                    "methods": ["igfold", "esmfold"],
                    "confidence_threshold": 0.7,
                    "batch_size": 2
                },
                "use_mock": self.use_mocks
            }
            
            # Save config to file
            config_path = evolution_dir / "evolution_config.json"
            with open(config_path, "w") as f:
                json.dump(config, f, indent=2)
            
            # Initialize hybrid evolutionary search
            start_time = time.time()
            hybrid_search = HybridEvolutionarySearch(str(config_path))
            
            # Run optimization
            logger.info(f"Starting optimization with {len(initial_antibodies)} initial antibodies...")
            optimized_antibodies = hybrid_search.optimize_antibody(
                initial_population=initial_antibodies,
                target_epitope=target_epitope
            )
            
            # Calculate total time
            total_time = time.time() - start_time
            
            # Record metrics
            self.metrics.add_metric("evolution_time", total_time)
            self.metrics.add_metric("num_optimized_antibodies", len(optimized_antibodies))
            
            if optimized_antibodies:
                best_ab = optimized_antibodies[0]  # Should be sorted by fitness
                self.metrics.add_metric("best_binding_energy", best_ab.binding_energy)
                self.metrics.add_metric("best_stability", best_ab.stability)
                self.metrics.add_metric("best_developability", best_ab.developability)
                
                logger.info(f"Best antibody: {best_ab.name}")
                logger.info(f"  - Binding energy: {best_ab.binding_energy:.2f}")
                logger.info(f"  - Stability: {best_ab.stability:.2f}")
                logger.info(f"  - Developability: {best_ab.developability:.2f}")
            
            # Create results dictionary
            results = {
                "success": True,
                "evolution_time": total_time,
                "generation_count": 3,
                "initial_population_size": len(initial_antibodies),
                "final_population_size": len(optimized_antibodies),
                "optimized_antibodies": [
                    {
                        "name": ab.name,
                        "binding_energy": getattr(ab, "binding_energy", 0),
                        "stability": getattr(ab, "stability", 0),
                        "developability": getattr(ab, "developability", 0),
                        "fitness": getattr(ab, "fitness", 0)
                    }
                    for ab in optimized_antibodies
                ]
            }
            
            # Save results
            with open(evolution_dir / "evolution_results.json", "w") as f:
                json.dump(results, f, indent=2)
                
            logger.info(f"Evolutionary search completed in {total_time:.2f}s")
            return results, optimized_antibodies
            
        except Exception as e:
            logger.error(f"Evolutionary search failed: {e}")
            return {"success": False, "error": str(e)}, []
    
    def test_cocktail_formulation(self, antibodies: List[Individual], toxins: Dict[str, Dict]) -> Dict[str, Any]:
        """
        Test the cocktail formulation strategy.
        
        Args:
            antibodies: List of antibodies to include in cocktail
            toxins: Dictionary of toxins to target
            
        Returns:
            Dictionary with test results
        """
        logger.info("Testing cocktail formulation...")
        
        # Create cocktail output directory
        cocktail_dir = self.output_dir / "cocktail_strategy"
        os.makedirs(cocktail_dir, exist_ok=True)
        
        try:
            # Prepare antibody data for cocktail optimizer
            formatted_antibodies = []
            for ab in antibodies:
                # Import required libraries
                import numpy as np
                import random
                
                # Create binding information (mock data)
                bindings = []
                for toxin_id, toxin in toxins.items():
                    # Create random but reasonable binding affinity
                    # Lower binding energy indicates stronger binding
                    binding_energy = getattr(ab, "binding_energy", -random.uniform(6, 12))
                    
                    # Convert binding energy to affinity (nM)
                    # Rough approximation: ΔG = -RT ln(1/Kd)
                    # So Kd = exp(-ΔG/RT)
                    # Using RT = 0.593 kcal/mol at room temperature
                    affinity = round(np.exp(-binding_energy/0.593), 2)
                    
                    bindings.append({
                        "toxin_id": toxin_id,
                        "affinity_nm": affinity,
                        "binding_energy": binding_energy
                    })
                
                # Numpy already imported above
                
                # Create formatted antibody with all required fields
                formatted_antibodies.append({
                    "id": ab.name,
                    "heavy_chain": ab.heavy_chain,
                    "light_chain": ab.light_chain,
                    "developability_score": getattr(ab, "developability", random.uniform(0.6, 0.9)),
                    "stability": getattr(ab, "stability", random.uniform(0.6, 0.9)),
                    "binds_to": bindings,
                    "is_broadly_neutralizing": random.choice([True, False]),
                    "targeting_strategy": "broadly_neutralizing" if random.random() > 0.6 else "specific"
                })
            
            # Initialize cocktail optimizer
            cocktail_optimizer = CocktailOptimizer()
            
            # Time the cocktail optimization
            start_time = time.time()
            
            # Design optimal cocktail
            cocktail = cocktail_optimizer.design_optimal_cocktail(formatted_antibodies, toxins)
            
            # Calculate total time
            total_time = time.time() - start_time
            
            # Record metrics
            self.metrics.add_metric("cocktail_optimization_time", total_time)
            self.metrics.add_metric("cocktail_antibody_count", len(cocktail.get("antibodies", [])))
            self.metrics.add_metric("cocktail_average_coverage", cocktail.get("average_coverage", 0))
            
            # Create results dictionary
            results = {
                "success": True,
                "optimization_time": total_time,
                "cocktail": cocktail
            }
            
            # Save results
            with open(cocktail_dir / "cocktail_results.json", "w") as f:
                json.dump(results, f, indent=2)
                
            logger.info(f"Cocktail formulation completed in {total_time:.2f}s")
            return results
            
        except Exception as e:
            logger.error(f"Cocktail formulation failed: {e}")
            return {"success": False, "error": str(e)}
    
    def visualize_results(self) -> None:
        """
        Generate visualizations of test results.
        """
        logger.info("Generating result visualizations...")
        
        # Create visualizations directory
        viz_dir = self.output_dir / "visualizations"
        os.makedirs(viz_dir, exist_ok=True)
        
        # Plot metrics
        self.metrics.plot_metrics(
            ["heavy_chain_prediction_time", "light_chain_prediction_time", "total_prediction_time"],
            output_path=str(viz_dir / "prediction_times.png"),
            title="Structure Prediction Times"
        )
        
        self.metrics.plot_metrics(
            ["heavy_chain_confidence", "light_chain_confidence"],
            output_path=str(viz_dir / "prediction_confidence.png"),
            title="Structure Prediction Confidence"
        )
        
        if "best_binding_energy" in self.metrics.metrics:
            self.metrics.plot_metrics(
                ["best_binding_energy", "best_stability", "best_developability"],
                output_path=str(viz_dir / "best_antibody_metrics.png"),
                title="Best Antibody Metrics"
            )
        
        # Export metrics summary
        metrics_json = self.metrics.export_to_json(str(viz_dir / "metrics_summary.json"))
        
        # Create summary report
        with open(viz_dir / "summary_report.md", "w") as f:
            f.write("# Phytovenomics Pipeline Demo Test Summary\n\n")
            f.write(f"Test completed at: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Structure Prediction\n")
            if "heavy_chain_prediction_time" in self.metrics.metrics:
                heavy_time = self.metrics.get_metric_average("heavy_chain_prediction_time")
                light_time = self.metrics.get_metric_average("light_chain_prediction_time")
                heavy_conf = self.metrics.get_metric_average("heavy_chain_confidence")
                light_conf = self.metrics.get_metric_average("light_chain_confidence")
                
                f.write(f"- Heavy chain prediction: {heavy_time:.2f}s, confidence: {heavy_conf:.2f}\n")
                f.write(f"- Light chain prediction: {light_time:.2f}s, confidence: {light_conf:.2f}\n")
            else:
                f.write("- Structure prediction was not tested or failed\n")
                
            f.write("\n## Evolutionary Optimization\n")
            if "evolution_time" in self.metrics.metrics:
                evol_time = self.metrics.get_metric_average("evolution_time")
                num_abs = self.metrics.get_metric_average("num_optimized_antibodies")
                f.write(f"- Optimization time: {evol_time:.2f}s\n")
                f.write(f"- Optimized antibodies: {num_abs:.0f}\n")
                
                if "best_binding_energy" in self.metrics.metrics:
                    best_energy = self.metrics.get_metric_average("best_binding_energy")
                    best_stab = self.metrics.get_metric_average("best_stability")
                    best_dev = self.metrics.get_metric_average("best_developability")
                    
                    f.write(f"- Best antibody metrics:\n")
                    f.write(f"  - Binding energy: {best_energy:.2f}\n")
                    f.write(f"  - Stability: {best_stab:.2f}\n")
                    f.write(f"  - Developability: {best_dev:.2f}\n")
            else:
                f.write("- Evolutionary optimization was not tested or failed\n")
                
            f.write("\n## Cocktail Formulation\n")
            if "cocktail_optimization_time" in self.metrics.metrics:
                cocktail_time = self.metrics.get_metric_average("cocktail_optimization_time")
                ab_count = self.metrics.get_metric_average("cocktail_antibody_count")
                coverage = self.metrics.get_metric_average("cocktail_average_coverage")
                
                f.write(f"- Cocktail formulation time: {cocktail_time:.2f}s\n")
                f.write(f"- Antibody count in cocktail: {ab_count:.0f}\n")
                f.write(f"- Average toxin coverage: {coverage:.2f}\n")
            else:
                f.write("- Cocktail formulation was not tested or failed\n")
        
        logger.info(f"Result visualizations created in {viz_dir}")
    
    def run_pipeline_test(self) -> Dict[str, Any]:
        """
        Run the complete pipeline test.
        
        Returns:
            Dictionary with test results
        """
        logger.info("Starting full pipeline demo test...")
        
        start_time = time.time()
        results = {"success": True}
        
        try:
            # 1. Create test data
            test_data = self.create_test_data()
            initial_antibodies = test_data["initial_antibodies"]
            target_epitope = test_data["epitope"]
            toxins = test_data["toxins"]
            
            # 2. Test structure prediction on first antibody
            structure_results = self.test_structure_prediction(initial_antibodies[0])
            results["structure_prediction"] = structure_results
            
            if not structure_results.get("success", False):
                logger.error("Structure prediction failed, but continuing with pipeline test")
            
            # 3. Test hybrid evolutionary search
            evolution_results, optimized_antibodies = self.test_hybrid_evolutionary_search(
                initial_antibodies, target_epitope
            )
            results["hybrid_evolution"] = evolution_results
            
            # If evolutionary search failed, use initial antibodies
            if not evolution_results.get("success", False) or not optimized_antibodies:
                logger.warning("Evolutionary search failed, using initial antibodies for cocktail formulation")
                optimized_antibodies = initial_antibodies
            
            # 4. Test cocktail formulation
            cocktail_results = self.test_cocktail_formulation(optimized_antibodies, toxins)
            results["cocktail_formulation"] = cocktail_results
            
            # 5. Generate visualizations
            self.visualize_results()
            
        except Exception as e:
            logger.error(f"Pipeline test failed: {e}")
            results["success"] = False
            results["error"] = str(e)
        
        # Record total time
        total_time = time.time() - start_time
        self.metrics.add_metric("total_pipeline_time", total_time)
        results["total_time"] = total_time
        
        # Save overall results
        with open(self.output_dir / "pipeline_test_results.json", "w") as f:
            json.dump(results, f, indent=2)
        
        if results["success"]:
            logger.info(f"Full pipeline test completed successfully in {total_time:.2f}s")
        else:
            logger.error(f"Pipeline test failed after {total_time:.2f}s")
        
        return results

def main():
    """Main function to run the demo test."""
    parser = argparse.ArgumentParser(description="Phytovenomics Pipeline Demo Test")
    parser.add_argument("--output-dir", default="demo_test_results", help="Output directory")
    parser.add_argument("--no-mocks", action="store_true", help="Disable mock implementations")
    parser.add_argument("--test", choices=["structure", "evolution", "cocktail", "all"], 
                       default="all", help="Specific test to run")
    parser.add_argument("--verbose", "-v", action="store_true", help="Increase verbosity")
    
    args = parser.parse_args()
    
    # Set verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Create demo test instance
    demo_test = PipelineDemoTest(
        output_dir=args.output_dir,
        use_mocks=not args.no_mocks
    )
    
    # Create test data
    test_data = demo_test.create_test_data()
    
    # Run specified test or all tests
    if args.test == "structure" or args.test == "all":
        demo_test.test_structure_prediction(test_data["initial_antibodies"][0])
        
    if args.test == "evolution" or args.test == "all":
        demo_test.test_hybrid_evolutionary_search(
            test_data["initial_antibodies"], test_data["epitope"]
        )
        
    if args.test == "cocktail" or args.test == "all":
        demo_test.test_cocktail_formulation(
            test_data["initial_antibodies"], test_data["toxins"]
        )
        
    if args.test == "all":
        demo_test.run_pipeline_test()
    
    # Generate visualizations
    demo_test.visualize_results()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())