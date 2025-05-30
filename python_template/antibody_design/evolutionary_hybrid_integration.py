#!/usr/bin/env python
# antibody_design/evolutionary_hybrid_integration.py

"""
Integration module for evolutionary search and hybrid structure prediction.
This module provides a unified pipeline that combines evolutionary optimization 
of antibody sequences with structure prediction to enhance the design process.
"""

import os
import sys
import time
import logging
import random
import json
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field

# Import modules from the project
from antibody_design.evolutionary_search import EvolutionarySearch, Individual, RosettaFoldClient
import utils.visualization as viz
from utils.validation_metrics import ValidationMetrics

# Try to import the structure prediction module
try:
    from src.structure_prediction.hybrid_prediction import HybridPredictionPipeline
except ImportError:
    logging.warning("Failed to import HybridPredictionPipeline, using mock implementation")
    
    # Create a mock class for testing if the real one is not available
    class HybridPredictionPipeline:
        def __init__(self, use_cache=True, cache_dir=None):
            self.use_cache = use_cache
            self.cache_dir = cache_dir
            
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

logger = logging.getLogger(__name__)


class HybridEvolutionarySearch:
    """
    Evolutionary search integrated with structure prediction for antibody design.
    This class extends the evolutionary search algorithm by incorporating structure prediction
    to evaluate the stability and developability of antibodies with higher accuracy.
    """
    
    def __init__(self, config_path_or_dict):
        """
        Initialize the hybrid evolutionary search algorithm.
        
        Args:
            config_path_or_dict: Path to configuration file or configuration dictionary
        """
        # Load configuration
        if isinstance(config_path_or_dict, str):
            with open(config_path_or_dict, 'r') as f:
                if config_path_or_dict.endswith('.json'):
                    self.config = json.load(f)
                elif config_path_or_dict.endswith(('.yaml', '.yml')):
                    import yaml
                    self.config = yaml.safe_load(f)
                else:
                    raise ValueError(f"Unsupported config file format: {config_path_or_dict}")
        else:
            self.config = config_path_or_dict
        
        # Extract general evolutionary search parameters
        self.population_size = self.config.get("population_size", 20)
        self.elite_size = self.config.get("elite_size", 2)
        self.max_generations = self.config.get("generations", 10)
        self.mutation_rate = self.config.get("mutation_rate", 0.1)
        self.crossover_rate = self.config.get("crossover_rate", 0.7)
        
        # Extract hybrid prediction specific parameters
        self.hybrid_config = self.config.get("hybrid_prediction", {})
        self.use_structure_prediction = self.hybrid_config.get("enabled", True)
        self.structure_methods = self.hybrid_config.get("methods", ["igfold", "esmfold"])
        self.confidence_threshold = self.hybrid_config.get("confidence_threshold", 0.7)
        self.batch_size = self.hybrid_config.get("batch_size", 4)
        
        # Set up output directory
        self.output_dir = self.config.get("output_dir", "hybrid_evolutionary_output")
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Set mock mode for testing
        self.use_mock = self.config.get("use_mock", os.getenv("USE_MOCK_EVALUATION", "True").lower() == "true")
        
        # Initialize the base evolutionary search
        self.evolutionary_search = EvolutionarySearch(self.config)
        
        # Initialize structure prediction pipeline
        if self.use_structure_prediction:
            # Initialize structure prediction pipeline
            # Note: Different parameters depending on what version of HybridPredictionPipeline is available
            try:
                self.structure_pipeline = HybridPredictionPipeline(
                    use_cache=True,
                    cache_dir=os.path.join(self.output_dir, "structure_cache"),
                    use_mock=self.use_mock,
                    models=self.structure_methods,
                    confidence_threshold=self.confidence_threshold
                )
            except TypeError:
                # Fallback for when use_mock is not a valid parameter
                logger.warning("Using HybridPredictionPipeline without mock parameter")
                self.structure_pipeline = HybridPredictionPipeline(
                    use_cache=True,
                    cache_dir=os.path.join(self.output_dir, "structure_cache")
                )
        
        # Initialize metrics tracker
        self.metrics = ValidationMetrics(self.output_dir)
        
        # Track the best individuals across generations
        self.best_individuals = []
        self._current_generation = 0
        
        logger.info(f"Initialized hybrid evolutionary search with population size {self.population_size}")
        logger.info(f"Using structure prediction: {self.use_structure_prediction}")
        logger.info(f"Using mock evaluations: {self.use_mock}")
    
    def _evaluate_individual_structure(self, individual: Individual) -> Individual:
        """
        Evaluate an individual using structure prediction to update stability and developability.
        
        Args:
            individual: Antibody individual to evaluate
            
        Returns:
            Evaluated individual with updated metrics
        """
        try:
            # Predict heavy chain structure
            logger.debug(f"Predicting structure for {individual.name} heavy chain")
            heavy_result = self.structure_pipeline.predict_antibody_structure(
                sequence=individual.heavy_chain,
                chain_type="heavy",
                use_ensemble=True
            )
            
            # Predict light chain structure
            logger.debug(f"Predicting structure for {individual.name} light chain")
            light_result = self.structure_pipeline.predict_antibody_structure(
                sequence=individual.light_chain,
                chain_type="light",
                use_ensemble=True
            )
            
            # Analyze structure quality
            heavy_quality = self.structure_pipeline.analyze_structure_quality(heavy_result)
            light_quality = self.structure_pipeline.analyze_structure_quality(light_result)
            
            # Update stability based on structure quality
            # Weight the heavy chain more than the light chain (60/40 split)
            heavy_stability = heavy_quality.get("confidence", 0) * 0.8
            light_stability = light_quality.get("confidence", 0) * 0.8
            
            # Combine into a single stability score
            combined_stability = 0.6 * heavy_stability + 0.4 * light_stability
            
            # Only update if the predicted stability is reliable
            if heavy_result.confidence > self.confidence_threshold and light_result.confidence > self.confidence_threshold:
                individual.stability = combined_stability
                
                # Also update developability score based on structural features
                # In a real implementation, we would analyze additional structural features
                individual.developability = 0.7 * combined_stability + 0.3 * random.uniform(0.6, 0.9)
            
            # Update individual fitness
            individual.update_fitness()
            
            # Track metrics
            self.metrics.add_metric("heavy_chain_confidence", heavy_result.confidence)
            self.metrics.add_metric("light_chain_confidence", light_result.confidence)
            self.metrics.add_metric("structure_stability", combined_stability)
            
            return individual
            
        except Exception as e:
            logger.error(f"Structure prediction failed for {individual.name}: {e}")
            # Keep original scores and return
            return individual
    
    def _evaluate_batch_with_structure(self, batch: List[Individual]) -> List[Individual]:
        """
        Evaluate a batch of individuals using structure prediction.
        
        Args:
            batch: List of individuals to evaluate
            
        Returns:
            List of evaluated individuals
        """
        evaluated_batch = []
        
        for individual in batch:
            try:
                # Evaluate using structure prediction
                evaluated_individual = self._evaluate_individual_structure(individual)
                evaluated_batch.append(evaluated_individual)
                
            except Exception as e:
                logger.error(f"Failed to evaluate {individual.name}: {e}")
                # Keep original individual
                evaluated_batch.append(individual)
        
        return evaluated_batch
    
    def _evaluate_population_with_structure(self, population: List[Individual]) -> List[Individual]:
        """
        Evaluate a population using structure prediction, processing in batches.
        
        Args:
            population: List of individuals to evaluate
            
        Returns:
            List of evaluated individuals
        """
        if not self.use_structure_prediction:
            logger.info("Structure prediction disabled, skipping structural evaluation")
            return population
        
        logger.info(f"Evaluating population with structure prediction (batch size: {self.batch_size})")
        
        # Process in batches to improve efficiency
        evaluated_population = []
        
        for i in range(0, len(population), self.batch_size):
            batch = population[i:i+self.batch_size]
            
            logger.info(f"Processing batch {i//self.batch_size + 1}/{(len(population)-1)//self.batch_size + 1} ({len(batch)} individuals)")
            
            # Evaluate batch
            evaluated_batch = self._evaluate_batch_with_structure(batch)
            evaluated_population.extend(evaluated_batch)
            
            logger.info(f"Batch {i//self.batch_size + 1} completed")
        
        return evaluated_population
    
    def optimize_antibody(self, initial_population: List[Individual], target_epitope: Dict) -> List[Individual]:
        """
        Optimize antibodies using evolutionary search with structure prediction integration.
        
        Args:
            initial_population: Initial population of antibody individuals
            target_epitope: Target epitope information
            
        Returns:
            List of optimized antibody individuals
        """
        logger.info(f"Starting hybrid evolutionary optimization with {len(initial_population)} initial antibodies")
        
        # Ensure we have enough initial individuals
        if len(initial_population) < self.population_size:
            logger.warning(f"Initial population size ({len(initial_population)}) is less than configured "
                         f"population size ({self.population_size}). Cloning to fill.")
            
            # Clone existing individuals to fill the population
            while len(initial_population) < self.population_size:
                to_clone = random.choice(initial_population)
                clone = Individual(
                    heavy_chain=to_clone.heavy_chain,
                    light_chain=to_clone.light_chain,
                    cdr_sequences=to_clone.cdr_sequences.copy(),
                    name=f"{to_clone.name}_clone{len(initial_population)}",
                    parent_ids=[to_clone.name],
                    generation=0
                )
                initial_population.append(clone)
        
        # First, evaluate the entire population with the base evolutionary search
        logger.info("Evaluating initial population with base metrics")
        population = self.evolutionary_search._evaluate_population(initial_population, target_epitope)
        
        # Then, evaluate with structure prediction
        if self.use_structure_prediction:
            logger.info("Enhancing evaluation with structure prediction")
            population = self._evaluate_population_with_structure(population)
        
        # Extract best initial individual
        best_initial = max(population, key=lambda x: x.fitness)
        self.best_individuals = [best_initial]
        
        logger.info(f"Initial best: {best_initial.name} with fitness {best_initial.fitness:.4f}")
        
        # Run evolutionary optimization for max_generations
        for gen in range(self.max_generations):
            self._current_generation = gen + 1
            logger.info(f"Starting generation {self._current_generation}/{self.max_generations}")
            
            # Create next generation with the base evolutionary search mechanism
            next_gen = self.evolutionary_search._create_next_generation(population, target_epitope)
            
            # Enhance evaluation with structure prediction
            if self.use_structure_prediction:
                # Only evaluate elite individuals and a random sample of the rest to save computational resources
                elite = sorted(next_gen, key=lambda x: x.fitness, reverse=True)[:self.elite_size]
                non_elite = sorted(next_gen, key=lambda x: x.fitness, reverse=True)[self.elite_size:]
                
                # Select candidates for structural evaluation (elite + sample of non-elite)
                struct_eval_count = min(len(next_gen), max(self.elite_size, self.batch_size))
                if len(non_elite) > 0 and struct_eval_count > self.elite_size:
                    sample_size = struct_eval_count - self.elite_size
                    non_elite_sample = random.sample(non_elite, min(sample_size, len(non_elite)))
                    struct_candidates = elite + non_elite_sample
                else:
                    struct_candidates = elite
                
                logger.info(f"Enhancing {len(struct_candidates)} individuals with structure prediction")
                enhanced_candidates = self._evaluate_population_with_structure(struct_candidates)
                
                # Replace the evaluated individuals in the population
                candidate_names = {ind.name for ind in struct_candidates}
                next_gen = [ind for ind in next_gen if ind.name not in candidate_names] + enhanced_candidates
            
            # Update population
            population = next_gen
            
            # Track best individual
            best_in_gen = max(population, key=lambda x: x.fitness)
            self.best_individuals.append(best_in_gen)
            
            # Track metrics
            avg_fitness = sum(ind.fitness for ind in population) / len(population)
            self.metrics.add_metric(f"generation_{gen+1}_best_fitness", best_in_gen.fitness)
            self.metrics.add_metric(f"generation_{gen+1}_avg_fitness", avg_fitness)
            self.metrics.add_metric(f"generation_{gen+1}_best_binding", best_in_gen.binding_energy)
            self.metrics.add_metric(f"generation_{gen+1}_best_stability", best_in_gen.stability)
            
            # Log progress
            logger.info(f"Generation {self._current_generation}: "
                      f"Best fitness = {best_in_gen.fitness:.4f}, "
                      f"Binding = {best_in_gen.binding_energy:.2f}, "
                      f"Stability = {best_in_gen.stability:.2f}, "
                      f"Avg fitness = {avg_fitness:.4f}")
            
            # Save generation data
            self._save_generation_data(population, gen + 1)
            
            # Create visualizations periodically
            if (gen + 1) % 2 == 0 or gen + 1 == self.max_generations:
                self._create_visualizations(population, gen + 1)
        
        # Sort final population by fitness
        final_population = sorted(population, key=lambda x: x.fitness, reverse=True)
        
        # Save optimization results
        self._save_optimization_results(final_population)
        
        logger.info(f"Hybrid optimization complete after {self.max_generations} generations")
        logger.info(f"Best individual: {final_population[0].name} with fitness {final_population[0].fitness:.4f}")
        
        return final_population
    
    def _save_generation_data(self, population: List[Individual], generation: int) -> None:
        """
        Save data for a generation.
        
        Args:
            population: Current population
            generation: Generation number
        """
        # Create directory for generation data
        gen_dir = os.path.join(self.output_dir, f"generation_{generation}")
        os.makedirs(gen_dir, exist_ok=True)
        
        # Save population data
        population_data = [
            {
                "name": ind.name,
                "heavy_chain": ind.heavy_chain,
                "light_chain": ind.light_chain,
                "binding_energy": ind.binding_energy,
                "stability": ind.stability,
                "developability": ind.developability,
                "manufacturability": getattr(ind, "manufacturability", 0.0),
                "fitness": ind.fitness,
                "parent_ids": ind.parent_ids,
                "generation": ind.generation
            }
            for ind in population
        ]
        
        # Write to file
        with open(os.path.join(gen_dir, "population.json"), "w") as f:
            json.dump(population_data, f, indent=2)
    
    def _save_optimization_results(self, final_population: List[Individual]) -> None:
        """
        Save the overall optimization results.
        
        Args:
            final_population: Final optimized population
        """
        # Create results directory
        results_dir = os.path.join(self.output_dir, "results")
        os.makedirs(results_dir, exist_ok=True)
        
        # Save best individuals across generations
        best_across_generations = [
            {
                "name": ind.name,
                "heavy_chain": ind.heavy_chain,
                "light_chain": ind.light_chain,
                "binding_energy": ind.binding_energy,
                "stability": ind.stability,
                "developability": ind.developability,
                "manufacturability": getattr(ind, "manufacturability", 0.0),
                "fitness": ind.fitness,
                "generation": ind.generation,
                "cdr_sequences": {k: v for k, v in ind.cdr_sequences.items()}
            }
            for ind in self.best_individuals
        ]
        
        with open(os.path.join(results_dir, "best_across_generations.json"), "w") as f:
            json.dump(best_across_generations, f, indent=2)
        
        # Save final population
        final_population_data = [
            {
                "name": ind.name,
                "heavy_chain": ind.heavy_chain,
                "light_chain": ind.light_chain,
                "binding_energy": ind.binding_energy,
                "stability": ind.stability,
                "developability": ind.developability,
                "manufacturability": getattr(ind, "manufacturability", 0.0),
                "fitness": ind.fitness,
                "generation": ind.generation,
                "cdr_sequences": {k: v for k, v in ind.cdr_sequences.items()}
            }
            for ind in final_population
        ]
        
        with open(os.path.join(results_dir, "final_population.json"), "w") as f:
            json.dump(final_population_data, f, indent=2)
        
        # Save metrics summary
        self.metrics.export_to_json(os.path.join(results_dir, "metrics.json"))
    
    def _create_visualizations(self, population: List[Individual], generation: int) -> None:
        """
        Create visualizations for a generation.
        
        Args:
            population: Current population
            generation: Generation number
        """
        # Create visualizations directory
        viz_dir = os.path.join(self.output_dir, "visualizations")
        os.makedirs(viz_dir, exist_ok=True)
        
        # Sort population by fitness
        sorted_pop = sorted(population, key=lambda x: x.fitness, reverse=True)
        
        try:
            # Create antibody visualizer
            antibody_viz = viz.AntibodyVisualizer(viz_dir)
            
            # Plot binding energy of top 10 individuals
            antibody_viz.plot_binding_energy(
                sorted_pop[:10],
                output_path=os.path.join(viz_dir, f"gen_{generation}_binding_energy.png")
            )
            
            # Plot metrics of top 5 individuals
            antibody_viz.plot_antibody_metrics(
                sorted_pop[:5],
                metrics=['binding_energy', 'stability', 'developability'],
                output_path=os.path.join(viz_dir, f"gen_{generation}_metrics.png")
            )
            
            # Plot evolution progress if we have best individuals tracked
            if self.best_individuals and generation > 1:
                # Extract generations and best fitness values
                gen_data = [{
                    'generation': ind.generation,
                    'best_fitness': ind.fitness,
                    'best_binding': ind.binding_energy,
                    'best_stability': ind.stability,
                    'best_developability': ind.developability
                } for ind in self.best_individuals]
                
                # Plot evolution progress
                antibody_viz.plot_evolution_progress(
                    gen_data,
                    metric='best_fitness',
                    output_path=os.path.join(viz_dir, f"evolution_fitness_{generation}.png")
                )
            
            # Plot CDR profile of the best antibody
            if sorted_pop:
                antibody_viz.plot_cdr_profile(
                    sorted_pop[0],
                    output_path=os.path.join(viz_dir, f"gen_{generation}_best_cdr_profile.png")
                )
            
            logger.info(f"Created visualizations for generation {generation}")
            
        except Exception as e:
            logger.error(f"Failed to create visualizations: {e}")
    
    def get_best_antibody(self) -> Optional[Individual]:
        """
        Get the best antibody from the optimization.
        
        Returns:
            Best antibody individual or None if optimization hasn't been run
        """
        if not self.best_individuals:
            logger.warning("No optimization has been run yet")
            return None
        
        return max(self.best_individuals, key=lambda x: x.fitness)


# Module test function
def test_hybrid_evolutionary_search():
    """
    Test function for the hybrid evolutionary search module.
    """
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    
    # Create config for testing
    config = {
        "population_size": 5,
        "generations": 3,
        "mutation_rate": 0.1,
        "crossover_rate": 0.7,
        "elite_size": 1,
        "output_dir": "test_hybrid_output",
        "hybrid_prediction": {
            "enabled": True,
            "methods": ["igfold", "esmfold"],
            "confidence_threshold": 0.7,
            "batch_size": 2
        },
        "use_mock": True
    }
    
    # Create sample antibodies
    from antibody_design.evolutionary_search import Individual
    
    initial_population = [
        Individual(
            heavy_chain="QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQRLEWMGWINAGNGNTKYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYCAR",
            light_chain="DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIKR",
            cdr_sequences={
                "H1": "GYTFTSYAMH", "H2": "WINAGNGNTKYSQKFQG", "H3": "VGLVRGAFDI",
                "L1": "RASQSISSYLN", "L2": "AASSLQS", "L3": "QQSYSTPPT"
            },
            name="Sample1"
        ),
        Individual(
            heavy_chain="EVQLVESGGGLVQPGGSLRLSCAASGFDFSRYDMSWVRQAPGKGLEWVSYISSSGSTIYYADSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAR",
            light_chain="EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSLTFGGGTKVEIK",
            cdr_sequences={
                "H1": "GFDFSRYDMS", "H2": "YISSSGSTIYYADSVKG", "H3": "DRGYSGSRGAFDI",
                "L1": "RASQSVSSSY", "L2": "GASSRAT", "L3": "QQYGSSLT"
            },
            name="Sample2"
        ),
    ]
    
    # Define target epitope
    target_epitope = {
        "id": "spike_rbd_001", 
        "sequence": "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVS",
        "structure": None,
        "organism": "SARS-CoV-2"
    }
    
    # Initialize hybrid evolutionary search
    search = HybridEvolutionarySearch(config)
    
    # Run optimization
    results = search.optimize_antibody(initial_population, target_epitope)
    
    # Print results
    logger.info("Optimization complete!")
    for i, antibody in enumerate(results[:3]):  # Show top 3
        logger.info(f"\nOptimal Antibody {i+1}:")
        logger.info(f"Name: {antibody.name}")
        logger.info(f"Binding Energy: {antibody.binding_energy:.2f}")
        logger.info(f"Stability: {antibody.stability:.2f}")
        logger.info(f"Developability: {antibody.developability:.2f}")
    
    return results


if __name__ == "__main__":
    test_hybrid_evolutionary_search()