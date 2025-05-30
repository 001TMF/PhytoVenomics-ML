#!/usr/bin/env python
# run_pipeline_tests.py

import os
import sys
import logging
import argparse
import time
import json
from datetime import datetime
import unittest
from pathlib import Path
import matplotlib.pyplot as plt
from typing import List, Dict, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("Phytovenomics.TestRunner")

# Import test module
from test_pipeline_integration import TestPipelineInteractions

class PipelineTestRunner:
    """
    Runner for executing and coordinating tests of the Phytovenomics ML platform pipeline.
    """
    
    def __init__(self, output_dir: str = "test_results"):
        """
        Initialize test runner.
        
        Args:
            output_dir: Directory for saving test results
        """
        self.output_dir = Path(output_dir)
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Create timestamp for this test run
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.results_dir = self.output_dir / f"test_run_{self.timestamp}"
        os.makedirs(self.results_dir, exist_ok=True)
        
        # Set up logging to file
        self.log_file = self.results_dir / "test_run.log"
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        logging.getLogger().addHandler(file_handler)
        
        # Test results
        self.test_results = []
        
    def run_test(self, test_name: str) -> Dict[str, Any]:
        """
        Run a specific test and record results.
        
        Args:
            test_name: Name of the test method to run
            
        Returns:
            Dictionary with test results
        """
        logger.info(f"Running test: {test_name}")
        
        # Create test suite with single test
        suite = unittest.TestSuite()
        suite.addTest(TestPipelineInteractions(test_name))
        
        # Run test and capture results
        start_time = time.time()
        test_runner = unittest.TextTestRunner(verbosity=2)
        test_result = test_runner.run(suite)
        end_time = time.time()
        
        # Process results
        success = test_result.wasSuccessful()
        errors = [str(err[1]) for err in test_result.errors]
        failures = [str(fail[1]) for fail in test_result.failures]
        skipped = [str(skip[1]) for skip in test_result.skipped]
        
        result_data = {
            "test_name": test_name,
            "success": success,
            "errors": errors,
            "failures": failures,
            "skipped": skipped,
            "runtime": round(end_time - start_time, 2)
        }
        
        # Add to results list
        self.test_results.append(result_data)
        
        # Log result
        if success:
            logger.info(f"Test {test_name} completed successfully in {result_data['runtime']} seconds")
        else:
            logger.error(f"Test {test_name} failed in {result_data['runtime']} seconds")
            if errors:
                logger.error(f"Errors: {errors}")
            if failures:
                logger.error(f"Failures: {failures}")
            if skipped:
                logger.warning(f"Skipped: {skipped}")
                
        return result_data
        
    def run_all_tests(self) -> None:
        """Run all pipeline integration tests in sequence."""
        logger.info("Starting complete pipeline integration test suite")
        
        # Define test sequence (in dependency order)
        test_sequence = [
            "test_toxin_database_loading",
            "test_epitope_discovery",
            "test_antibody_generator",
            "test_affinity_optimizer",
            "test_evolutionary_search",
            "test_hybrid_prediction_pipeline",
            "test_hybrid_evolutionary_search",
            "test_cocktail_formulator",
            "test_validation_metrics",
            "test_end_to_end_pipeline"
        ]
        
        # Run each test
        for test_name in test_sequence:
            self.run_test(test_name)
            
        # Generate test report
        self.generate_report()
        
    def run_component_tests(self, components: List[str]) -> None:
        """
        Run tests for specific components.
        
        Args:
            components: List of component names to test
        """
        # Map component names to test methods
        component_map = {
            "toxin_database": "test_toxin_database_loading",
            "epitope_discovery": "test_epitope_discovery",
            "antibody_generator": "test_antibody_generator",
            "affinity_optimizer": "test_affinity_optimizer",
            "evolutionary_search": "test_evolutionary_search",
            "hybrid_prediction": "test_hybrid_prediction_pipeline",
            "hybrid_evolutionary": "test_hybrid_evolutionary_search",
            "cocktail_formulator": "test_cocktail_formulator",
            "validation_metrics": "test_validation_metrics",
            "pipeline": "test_end_to_end_pipeline"
        }
        
        # Validate component names
        unknown_components = [c for c in components if c not in component_map]
        if unknown_components:
            logger.warning(f"Unknown components: {unknown_components}")
            
        # Get tests to run
        tests_to_run = [component_map[c] for c in components if c in component_map]
        
        if not tests_to_run:
            logger.error("No valid components specified for testing")
            return
            
        # Run specified tests
        for test_name in tests_to_run:
            self.run_test(test_name)
            
        # Generate test report
        self.generate_report()
        
    def generate_report(self) -> None:
        """Generate test report and visualization."""
        # Save results as JSON
        results_file = self.results_dir / "test_results.json"
        with open(results_file, 'w') as f:
            json.dump(self.test_results, f, indent=2)
            
        # Generate metrics
        successful = sum(1 for r in self.test_results if r["success"])
        failed = len(self.test_results) - successful
        
        # Create summary report
        summary_file = self.results_dir / "test_summary.md"
        with open(summary_file, 'w') as f:
            f.write("# Phytovenomics Pipeline Integration Test Summary\n\n")
            f.write(f"Test Run: {self.timestamp}\n\n")
            f.write(f"**Overall Results**: {successful}/{len(self.test_results)} tests passed\n\n")
            
            f.write("## Test Results\n\n")
            f.write("| Test | Status | Runtime (s) | Issues |\n")
            f.write("|------|--------|------------|--------|\n")
            
            for result in self.test_results:
                status = "✅ PASS" if result["success"] else "❌ FAIL"
                issues = []
                if result["errors"]:
                    issues.extend([f"Error: {e}" for e in result["errors"]])
                if result["failures"]:
                    issues.extend([f"Failure: {f}" for f in result["failures"]])
                if result["skipped"]:
                    issues.extend([f"Skipped: {s}" for s in result["skipped"]])
                
                issue_text = "<br>".join(issues) if issues else "None"
                f.write(f"| {result['test_name']} | {status} | {result['runtime']} | {issue_text} |\n")
                
            # Add overall metrics
            f.write("\n## Performance Metrics\n\n")
            total_time = sum(r["runtime"] for r in self.test_results)
            f.write(f"Total runtime: {total_time:.2f} seconds\n\n")
            f.write(f"Average test time: {total_time/len(self.test_results):.2f} seconds\n\n")
            
            # Add log file location
            f.write("\n## Test Artifacts\n\n")
            f.write(f"- Log file: {self.log_file}\n")
            f.write(f"- Results JSON: {results_file}\n")
        
        # Create visualization
        self._create_test_visualization()
        
        logger.info(f"Test report generated at {self.results_dir}")
        logger.info(f"Summary: {successful}/{len(self.test_results)} tests passed")
        
    def _create_test_visualization(self) -> None:
        """Create visualization of test results."""
        # Extract data
        test_names = [r["test_name"].replace("test_", "").replace("_", " ").title() for r in self.test_results]
        runtimes = [r["runtime"] for r in self.test_results]
        success = [r["success"] for r in self.test_results]
        
        # Create bar chart of test runtimes, colored by success/failure
        plt.figure(figsize=(12, 8))
        bar_colors = ['green' if s else 'red' for s in success]
        
        # Plot horizontal bars
        bars = plt.barh(test_names, runtimes, color=bar_colors)
        plt.xlabel('Runtime (seconds)')
        plt.title('Phytovenomics Pipeline Test Results')
        plt.grid(axis='x', linestyle='--', alpha=0.7)
        
        # Add runtime labels to bars
        for bar, runtime in zip(bars, runtimes):
            plt.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2, 
                     f'{runtime:.2f}s', va='center')
        
        # Add success/failure counts in legend
        successful = sum(success)
        plt.legend([
            f'Passed ({successful}/{len(success)})', 
            f'Failed ({len(success)-successful}/{len(success)})'
        ], loc='best')
        
        # Save figure
        plt.tight_layout()
        plt.savefig(self.results_dir / 'test_results.png', dpi=300)
        plt.close()

def main():
    """Main function to run tests."""
    parser = argparse.ArgumentParser(description="Phytovenomics Pipeline Test Runner")
    
    # Test selection options
    test_group = parser.add_mutually_exclusive_group(required=True)
    test_group.add_argument("--all", action="store_true", help="Run all tests")
    test_group.add_argument("--components", nargs="+", help="Specific components to test")
    
    # Output directory
    parser.add_argument("--output-dir", default="test_results", help="Output directory for test results")
    
    # Verbosity
    parser.add_argument("--verbose", "-v", action="count", default=0, help="Increase verbosity")
    
    args = parser.parse_args()
    
    # Set verbosity
    if args.verbose >= 2:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.verbose == 1:
        logging.getLogger().setLevel(logging.INFO)
    
    # Create test runner
    runner = PipelineTestRunner(output_dir=args.output_dir)
    
    # Run tests
    if args.all:
        runner.run_all_tests()
    elif args.components:
        runner.run_component_tests(args.components)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())