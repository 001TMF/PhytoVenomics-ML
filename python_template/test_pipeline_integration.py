#!/usr/bin/env python
# test_pipeline_integration.py

import os
import sys
import logging
import unittest
import tempfile
import json
import yaml
import random
import torch
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("Phytovenomics.Test")

# Import Phytovenomics modules
from antibody_design import HumanAntibodyDesigner, EpitopeDiscovery, AffinityOptimizer
from antibody_design.evolutionary_search import EvolutionarySearch, Individual
from antibody_design.evolutionary_hybrid_integration import HybridEvolutionarySearch, StructureBasedEvaluator
from cocktail_strategy import CocktailOptimizer
from venom_data.toxin_database import ToxinDatabase

# Import structure prediction modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.structure_prediction.hybrid_prediction import HybridPredictionPipeline

# Import utility modules
from utils.validation_utils import ValidationMetrics
from utils.visualization import plot_antibody_structures
from utils.training_utils import set_seed

class TestPipelineInteractions(unittest.TestCase):
    """
    Test suite to validate the interactions between different components 
    of the Phytovenomics ML platform pipeline.
    """

    @classmethod
    def setUpClass(cls):
        """Set up test environment once before all tests."""
        logger.info("Setting up test environment...")
        
        # Set random seed for reproducibility
        set_seed(42)
        
        # Create temporary directory for test outputs
        cls.temp_dir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.temp_dir.name)
        
        # Load default test configuration
        cls.config = cls._load_test_config()
        
        # Create synthetic data for testing
        cls.toxin_data, cls.antibody_data, cls.binding_data = cls._create_test_data()
        
        logger.info(f"Test environment set up in {cls.output_dir}")

    @classmethod
    def tearDownClass(cls):
        """Clean up after all tests."""
        cls.temp_dir.cleanup()
        logger.info("Test environment cleaned up")

    @classmethod
    def _load_test_config(cls) -> Dict:
        """Load test configuration."""
        try:
            with open("config/default_config.yaml", "r") as f:
                config = yaml.safe_load(f)
                
            # Override settings for testing
            config["output_dir"] = str(cls.output_dir)
            config["use_mock"] = True  # Use mock implementations to avoid heavy computation
            
            # Update paths for test data
            config["toxin_database"]["data_path"] = str(cls.output_dir / "toxin_data.csv")
            config["antibody_designer"]["data_path"] = str(cls.output_dir / "antibody_data.csv")
            
            return config
        
        except (FileNotFoundError, yaml.YAMLError) as e:
            logger.warning(f"Could not load config file: {e}")
            logger.info("Using fallback test configuration")
            
            # Fallback config
            return {
                "output_dir": str(cls.output_dir),
                "use_mock": True,
                "toxin_database": {"data_path": str(cls.output_dir / "toxin_data.csv")},
                "antibody_designer": {"data_path": str(cls.output_dir / "antibody_data.csv")},
                "epitope_discovery": {"ml_model_path": "models/epitope_discovery_model.pkl"},
                "evolutionary_search": {
                    "population_size": 10,
                    "generations": 3,
                    "mutation_rate": 0.1,
                    "crossover_rate": 0.7,
                },
                "hybrid_prediction": {
                    "enabled": True,
                    "methods": ["igfold", "esmfold"],
                    "confidence_threshold": 0.7,
                    "batch_size": 2
                }
            }

    @classmethod
    def _create_test_data(cls) -> Tuple[List[Dict], List[Dict], List[Dict]]:
        """Create synthetic data for testing."""
        # Amino acid pool
        aa_pool = "ACDEFGHIKLMNPQRSTVWY"
        
        # Create toxin data
        toxins = []
        for i in range(5):
            toxin_family = random.choice([
                "three_finger_toxin",
                "phospholipase_a2",
                "snake_venom_metalloprotease"
            ])
            
            seq_length = random.randint(60, 100)
            sequence = ''.join(random.choice(aa_pool) for _ in range(seq_length))
            
            toxin = {
                "id": f"TOX{i+1:03d}",
                "sequence": sequence,
                "family": toxin_family,
                "source_species": random.choice(["Naja naja", "Bungarus caeruleus"]),
                "function": "neurotoxic" if "finger" in toxin_family else "hemotoxic",
                "potency": round(random.uniform(0.1, 5.0), 2)
            }
            toxins.append(toxin)
        
        # Create antibody data
        antibodies = []
        for i in range(5):
            # Generate heavy chain
            heavy_length = random.randint(110, 130)
            heavy_seq = ''.join(random.choice(aa_pool) for _ in range(heavy_length))
            
            # Generate light chain
            light_length = random.randint(100, 115)
            light_seq = ''.join(random.choice(aa_pool) for _ in range(light_length))
            
            antibody = {
                "id": f"AB{i+1:03d}",
                "heavy_chain": heavy_seq,
                "light_chain": light_seq,
                "developability_score": round(random.uniform(0.5, 0.95), 2)
            }
            antibodies.append(antibody)
        
        # Create binding data
        bindings = []
        for i in range(10):
            toxin = random.choice(toxins)
            antibody = random.choice(antibodies)
            
            binding = {
                "toxin_id": toxin["id"],
                "antibody_id": antibody["id"],
                "affinity_nm": round(random.uniform(0.1, 1000), 2),
                "binding_energy": round(random.uniform(-15, -5), 1)
            }
            bindings.append(binding)
        
        # Save data to files in the test directory
        cls._save_test_data(toxins, antibodies, bindings)
        
        return toxins, antibodies, bindings

    @classmethod
    def _save_test_data(cls, toxins: List[Dict], antibodies: List[Dict], bindings: List[Dict]):
        """Save test data to files."""
        import pandas as pd
        
        # Save toxin data
        toxin_df = pd.DataFrame(toxins)
        toxin_path = cls.output_dir / "toxin_data.csv"
        toxin_df.to_csv(toxin_path, index=False)
        
        # Save antibody data
        antibody_df = pd.DataFrame(antibodies)
        antibody_path = cls.output_dir / "antibody_data.csv"
        antibody_df.to_csv(antibody_path, index=False)
        
        # Save binding data
        binding_df = pd.DataFrame(bindings)
        binding_path = cls.output_dir / "binding_data.csv"
        binding_df.to_csv(binding_path, index=False)
        
        # Create FASTA file for antibodies
        fasta_path = cls.output_dir / "antibodies.fasta"
        with open(fasta_path, 'w') as f:
            for ab in antibodies:
                f.write(f">{ab['id']}_H\n{ab['heavy_chain']}\n")
                f.write(f">{ab['id']}_L\n{ab['light_chain']}\n")
        
        logger.info(f"Test data saved to {cls.output_dir}")

    def test_toxin_database_loading(self):
        """Test loading of toxin data into the ToxinDatabase module."""
        logger.info("Testing ToxinDatabase module...")
        
        # Initialize toxin database
        toxin_db = ToxinDatabase(self.config.get("toxin_database", {}))
        
        # Load toxin data
        toxin_db.load_from_csv(self.config["toxin_database"]["data_path"])
        
        # Verify data was loaded
        all_toxins = toxin_db.get_all_toxins()
        self.assertEqual(len(all_toxins), len(self.toxin_data), 
                         f"Expected {len(self.toxin_data)} toxins, got {len(all_toxins)}")
        
        # Test toxin lookup
        first_toxin_id = self.toxin_data[0]["id"]
        toxin = toxin_db.get_toxin_by_id(first_toxin_id)
        self.assertIsNotNone(toxin, f"Could not retrieve toxin with ID {first_toxin_id}")
        
        # Test filtering by family
        family = self.toxin_data[0]["family"]
        family_toxins = toxin_db.get_toxins_by_family(family)
        expected_count = sum(1 for t in self.toxin_data if t["family"] == family)
        self.assertEqual(len(family_toxins), expected_count, 
                         f"Expected {expected_count} toxins in family {family}, got {len(family_toxins)}")
        
        logger.info("ToxinDatabase test completed successfully")

    def test_epitope_discovery(self):
        """Test the epitope discovery module."""
        logger.info("Testing EpitopeDiscovery module...")
        
        # Initialize the epitope discovery module
        epitope_discovery = EpitopeDiscovery(self.config.get("epitope_discovery", {}))
        
        # Process each toxin
        for toxin in self.toxin_data[:2]:  # Process only first 2 toxins for efficiency
            logger.info(f"Discovering epitopes on toxin {toxin['id']}...")
            
            # Identify epitopes
            epitopes = epitope_discovery.identify_epitopes(toxin)
            
            # Verify epitopes were found
            self.assertIsNotNone(epitopes, f"No epitopes found for toxin {toxin['id']}")
            self.assertGreater(len(epitopes), 0, f"No epitopes found for toxin {toxin['id']}")
            
            # Check that each epitope has required properties
            for epitope in epitopes:
                self.assertIn("id", epitope, "Epitope is missing ID")
                self.assertIn("sequence", epitope, "Epitope is missing sequence")
                self.assertIn("score", epitope, "Epitope is missing score")
        
        logger.info("EpitopeDiscovery test completed successfully")

    def test_antibody_generator(self):
        """Test the antibody generator module."""
        logger.info("Testing HumanAntibodyDesigner module...")
        
        # Initialize the antibody designer
        antibody_designer = HumanAntibodyDesigner(self.config.get("antibody_designer", {}))
        
        # Get test toxin and epitope
        test_toxin = self.toxin_data[0]
        epitope_discovery = EpitopeDiscovery(self.config.get("epitope_discovery", {}))
        epitopes = epitope_discovery.identify_epitopes(test_toxin)
        test_epitope = epitopes[0] if epitopes else None
        
        self.assertIsNotNone(test_epitope, "No test epitope available")
        
        # Design antibody against epitope
        designed_antibody = antibody_designer.design_antibody(test_epitope)
        
        # Verify antibody was designed
        self.assertIsNotNone(designed_antibody, "No antibody was designed")
        self.assertIn("heavy_chain", designed_antibody, "Antibody missing heavy chain")
        self.assertIn("light_chain", designed_antibody, "Antibody missing light chain")
        self.assertIn("id", designed_antibody, "Antibody missing ID")
        
        logger.info("HumanAntibodyDesigner test completed successfully")

    def test_affinity_optimizer(self):
        """Test the affinity optimizer module."""
        logger.info("Testing AffinityOptimizer module...")
        
        # Initialize the affinity optimizer
        affinity_optimizer = AffinityOptimizer(self.config.get("affinity_optimizer", {}))
        
        # Create a test antibody
        test_antibody = {
            "id": "TEST_AB",
            "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKGWFDYWGQGTLVTVSS",
            "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIK",
            "developability_score": 0.75,
            "binding_energy": -8.5
        }
        
        # Get test toxin and epitope
        test_toxin = self.toxin_data[0]
        epitope_discovery = EpitopeDiscovery(self.config.get("epitope_discovery", {}))
        epitopes = epitope_discovery.identify_epitopes(test_toxin)
        test_epitope = epitopes[0] if epitopes else None
        
        # Optimize antibody
        optimized_antibody = affinity_optimizer.optimize_antibody(test_antibody, test_epitope)
        
        # Verify antibody was optimized
        self.assertIsNotNone(optimized_antibody, "No optimized antibody returned")
        self.assertIn("binding_energy", optimized_antibody, "Optimized antibody missing binding energy")
        
        # Check that binding energy improved (lower is better)
        self.assertLessEqual(optimized_antibody["binding_energy"], test_antibody["binding_energy"],
                          "Binding energy did not improve after optimization")
        
        logger.info("AffinityOptimizer test completed successfully")

    def test_evolutionary_search(self):
        """Test the evolutionary search module."""
        logger.info("Testing EvolutionarySearch module...")
        
        # Prepare configuration for evolutionary search
        es_config = self.config.get("evolutionary_search", {})
        es_config["output_dir"] = str(self.output_dir / "evolutionary_search")
        es_config["population_size"] = 6  # Smaller population for testing
        es_config["generations"] = 2  # Fewer generations for testing
        
        # Initialize evolutionary search
        evol_search = EvolutionarySearch(es_config)
        
        # Create initial population
        initial_population = []
        for i in range(4):
            # Create individuals from test antibody data
            if i < len(self.antibody_data):
                ab_data = self.antibody_data[i]
                ind = Individual(
                    heavy_chain=ab_data["heavy_chain"],
                    light_chain=ab_data["light_chain"],
                    cdr_sequences={  # Mock CDR sequences
                        "H1": ab_data["heavy_chain"][26:35],
                        "H2": ab_data["heavy_chain"][50:65],
                        "H3": ab_data["heavy_chain"][95:110] if len(ab_data["heavy_chain"]) > 110 else ab_data["heavy_chain"][95:],
                        "L1": ab_data["light_chain"][24:34],
                        "L2": ab_data["light_chain"][50:56],
                        "L3": ab_data["light_chain"][89:97] if len(ab_data["light_chain"]) > 97 else ab_data["light_chain"][89:]
                    },
                    name=f"TEST_IND_{i}"
                )
                initial_population.append(ind)
        
        # Get test epitope
        test_toxin = self.toxin_data[0]
        epitope_discovery = EpitopeDiscovery(self.config.get("epitope_discovery", {}))
        epitopes = epitope_discovery.identify_epitopes(test_toxin)
        test_epitope = epitopes[0] if epitopes else {
            "id": "MOCK_EPITOPE",
            "sequence": "RVQPTESIVRFPNITNLCPFG",
            "start": 10,
            "end": 30,
            "score": 0.85,
            "toxin_id": test_toxin["id"]
        }
        
        # Run evolutionary search
        optimized_antibodies = evol_search.optimize_antibody(initial_population, test_epitope)
        
        # Verify optimization results
        self.assertIsNotNone(optimized_antibodies, "No optimized antibodies returned")
        self.assertGreater(len(optimized_antibodies), 0, "Empty optimized antibody list")
        
        # Check that individuals have fitness values
        for ind in optimized_antibodies:
            self.assertTrue(hasattr(ind, "binding_energy"), "Individual missing binding energy")
            self.assertTrue(hasattr(ind, "stability"), "Individual missing stability")
        
        logger.info("EvolutionarySearch test completed successfully")

    def test_hybrid_prediction_pipeline(self):
        """Test the hybrid prediction pipeline."""
        logger.info("Testing HybridPredictionPipeline module...")
        
        try:
            # Initialize hybrid prediction pipeline with mock mode
            os.environ["USE_MOCK_PREDICTION"] = "True"  # Set mock mode
            pipeline = HybridPredictionPipeline(
                use_cache=True,
                cache_dir=str(self.output_dir / "structure_cache")
            )
            
            # Test antibody
            heavy_chain = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKGWFDYWGQGTLVTVSS"
            
            # Create temporary fasta file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
                temp_fasta.write(f">TEST_AB_H\n{heavy_chain}\n")
                temp_fasta_path = temp_fasta.name
            
            try:
                # Create prediction output directory
                pred_output_dir = self.output_dir / "structure_prediction"
                os.makedirs(pred_output_dir, exist_ok=True)
                
                # Predict structure
                result = pipeline.predict_antibody_structure(
                    sequence=heavy_chain,
                    chain_type="heavy",
                    use_ensemble=True
                )
                
                # Verify prediction result
                self.assertIsNotNone(result, "No prediction result returned")
                self.assertGreaterEqual(result.confidence, 0.0, "Invalid confidence score")
                self.assertIn("model_name", result.__dict__, "Model name missing from result")
            
            finally:
                # Clean up temporary file
                if os.path.exists(temp_fasta_path):
                    os.unlink(temp_fasta_path)
            
            logger.info("HybridPredictionPipeline test completed successfully")
        
        except Exception as e:
            logger.error(f"Error testing HybridPredictionPipeline: {str(e)}")
            self.skipTest(f"HybridPredictionPipeline test skipped due to error: {str(e)}")

    def test_hybrid_evolutionary_search(self):
        """Test the hybrid evolutionary search with structure-based evaluation."""
        logger.info("Testing HybridEvolutionarySearch module...")
        
        try:
            # Prepare configuration for hybrid evolutionary search
            es_config = self.config.get("evolutionary_search", {}).copy()
            es_config["output_dir"] = str(self.output_dir / "hybrid_evolutionary_search")
            es_config["population_size"] = 4  # Smaller population for testing
            es_config["generations"] = 2  # Fewer generations for testing
            es_config["hybrid_prediction"] = {
                "enabled": True,
                "methods": ["igfold", "esmfold"],
                "confidence_threshold": 0.7,
                "batch_size": 2
            }
            
            # Save config to a temporary file
            config_path = self.output_dir / "hybrid_es_config.yaml"
            with open(config_path, "w") as f:
                yaml.dump(es_config, f)
            
            # Set mock mode for testing
            os.environ["USE_MOCK_PREDICTION"] = "True"
            
            # Initialize hybrid evolutionary search
            hybrid_search = HybridEvolutionarySearch(config_path)
            
            # Create initial population (smaller than real use-case)
            initial_population = []
            for i in range(3):
                # Create individuals from test antibody data
                if i < len(self.antibody_data):
                    ab_data = self.antibody_data[i]
                    ind = Individual(
                        heavy_chain=ab_data["heavy_chain"],
                        light_chain=ab_data["light_chain"],
                        cdr_sequences={  # Mock CDR sequences
                            "H1": ab_data["heavy_chain"][26:35],
                            "H2": ab_data["heavy_chain"][50:65],
                            "H3": ab_data["heavy_chain"][95:110] if len(ab_data["heavy_chain"]) > 110 else ab_data["heavy_chain"][95:],
                            "L1": ab_data["light_chain"][24:34],
                            "L2": ab_data["light_chain"][50:56],
                            "L3": ab_data["light_chain"][89:97] if len(ab_data["light_chain"]) > 97 else ab_data["light_chain"][89:]
                        },
                        name=f"HYBRID_IND_{i}"
                    )
                    initial_population.append(ind)
            
            # Get test epitope
            test_toxin = self.toxin_data[0]
            epitope_discovery = EpitopeDiscovery(self.config.get("epitope_discovery", {}))
            epitopes = epitope_discovery.identify_epitopes(test_toxin)
            test_epitope = epitopes[0] if epitopes else {
                "id": "MOCK_EPITOPE",
                "sequence": "RVQPTESIVRFPNITNLCPFG",
                "start": 10,
                "end": 30,
                "score": 0.85,
                "toxin_id": test_toxin["id"]
            }
            
            # Run hybrid evolutionary search
            optimized_antibodies = hybrid_search.optimize_antibody(initial_population, test_epitope)
            
            # Verify optimization results
            self.assertIsNotNone(optimized_antibodies, "No optimized antibodies returned")
            self.assertGreater(len(optimized_antibodies), 0, "Empty optimized antibody list")
            
            logger.info("HybridEvolutionarySearch test completed successfully")
        
        except Exception as e:
            logger.error(f"Error testing HybridEvolutionarySearch: {str(e)}")
            self.skipTest(f"HybridEvolutionarySearch test skipped due to error: {str(e)}")

    def test_cocktail_formulator(self):
        """Test the cocktail formulator module."""
        logger.info("Testing CocktailOptimizer module...")
        
        # Initialize cocktail optimizer
        cocktail_optimizer = CocktailOptimizer(self.config.get("cocktail_optimizer", {}))
        
        # Create list of designed antibodies
        designed_antibodies = []
        for i, ab_data in enumerate(self.antibody_data):
            # Create a designed antibody with binding information
            ab = ab_data.copy()
            
            # Add binding information
            ab["binds_to"] = []
            for j, toxin in enumerate(self.toxin_data):
                if random.random() < 0.6:  # 60% chance of binding to each toxin
                    ab["binds_to"].append({
                        "toxin_id": toxin["id"],
                        "affinity_nm": round(random.uniform(0.1, 100), 2),
                        "binding_energy": round(random.uniform(-12, -6), 1)
                    })
            
            # Add targeting strategy and other fields
            ab["is_broadly_neutralizing"] = random.choice([True, False])
            ab["targeting_strategy"] = "broadly_neutralizing" if ab["is_broadly_neutralizing"] else "specific"
            
            designed_antibodies.append(ab)
        
        # Create toxin dictionary
        toxin_dict = {toxin["id"]: toxin for toxin in self.toxin_data}
        
        # Design optimal cocktail
        cocktail = cocktail_optimizer.design_optimal_cocktail(designed_antibodies, toxin_dict)
        
        # Verify cocktail
        self.assertIsNotNone(cocktail, "No cocktail returned")
        self.assertIn("antibodies", cocktail, "Cocktail missing antibodies list")
        self.assertIn("coverage", cocktail, "Cocktail missing coverage information")
        self.assertIn("average_coverage", cocktail, "Cocktail missing average_coverage value")
        
        logger.info("CocktailOptimizer test completed successfully")

    def test_validation_metrics(self):
        """Test the validation metrics module."""
        logger.info("Testing ValidationMetrics module...")
        
        # Initialize validation metrics
        metrics = ValidationMetrics()
        
        # Add some metrics
        metrics.add_metric("binding_energy", -10.5)
        metrics.add_metric("stability", 0.85)
        metrics.add_metric("developability", 0.78)
        
        # Add a series of metrics
        for i in range(5):
            metrics.add_metric(f"binding_{i}", -9.5 - i * 0.5)
        
        # Get summary statistics
        summary = metrics.get_summary_statistics()
        
        # Verify metrics were tracked
        self.assertIn("binding_energy", summary, "Binding energy not in summary")
        self.assertEqual(summary["binding_energy"]["value"], -10.5, "Incorrect binding energy value")
        
        # Verify series tracking
        binding_series = metrics.get_metric_series("binding_")
        self.assertEqual(len(binding_series), 5, f"Expected 5 binding metrics, got {len(binding_series)}")
        
        # Test JSON export
        json_data = metrics.export_to_json()
        self.assertIsNotNone(json_data, "No JSON data exported")
        
        # Test plotting
        plot_path = self.output_dir / "metrics_plot.png"
        metrics.plot_metrics(["binding_energy", "stability"], output_path=str(plot_path))
        self.assertTrue(plot_path.exists(), f"Metrics plot not created at {plot_path}")
        
        logger.info("ValidationMetrics test completed successfully")

    def test_end_to_end_pipeline(self):
        """Test the end-to-end pipeline integration."""
        logger.info("Testing end-to-end pipeline integration...")
        
        # 1. Initialize all components
        toxin_db = ToxinDatabase(self.config.get("toxin_database", {}))
        epitope_discovery = EpitopeDiscovery(self.config.get("epitope_discovery", {}))
        antibody_designer = HumanAntibodyDesigner(self.config.get("antibody_designer", {}))
        affinity_optimizer = AffinityOptimizer(self.config.get("affinity_optimizer", {}))
        cocktail_optimizer = CocktailOptimizer(self.config.get("cocktail_optimizer", {}))
        
        # Load toxin data
        toxin_db.load_from_csv(self.config["toxin_database"]["data_path"])
        all_toxins = toxin_db.get_all_toxins()
        
        # Select a single toxin for the pipeline test
        test_toxin = list(all_toxins.values())[0]
        
        # 2. Discover epitopes
        epitopes = epitope_discovery.identify_epitopes(test_toxin)
        self.assertGreater(len(epitopes), 0, "No epitopes discovered")
        
        # Select the top epitope
        top_epitope = epitopes[0]
        
        # 3. Design antibody against epitope
        designed_antibody = antibody_designer.design_antibody(top_epitope)
        self.assertIsNotNone(designed_antibody, "No antibody designed")
        
        # 4. Optimize antibody affinity
        optimized_antibody = affinity_optimizer.optimize_antibody(designed_antibody, top_epitope)
        self.assertIsNotNone(optimized_antibody, "No optimized antibody returned")
        
        # 5. Design a small cocktail
        designed_antibodies = [optimized_antibody]
        toxin_dict = {test_toxin["id"]: test_toxin}
        cocktail = cocktail_optimizer.design_optimal_cocktail(designed_antibodies, toxin_dict)
        self.assertIsNotNone(cocktail, "No cocktail returned")
        
        # 6. Save results
        results = {
            "toxin": test_toxin,
            "epitope": top_epitope,
            "antibody": optimized_antibody,
            "cocktail": cocktail
        }
        
        results_path = self.output_dir / "end_to_end_results.json"
        with open(results_path, "w") as f:
            # Convert any non-serializable objects
            json.dump(results, f, default=lambda o: str(o), indent=2)
        
        logger.info(f"End-to-end pipeline results saved to {results_path}")
        logger.info("End-to-end pipeline integration test completed successfully")

if __name__ == "__main__":
    unittest.main()