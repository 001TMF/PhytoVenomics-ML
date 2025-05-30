#!/usr/bin/env python
# utils/data_utils.py

import logging
import json
import os
from typing import Dict, List, Any, Optional, Union
import csv
import pandas as pd
import numpy as np
from pathlib import Path

logger = logging.getLogger("Phytovenomics.DataUtils")

class DataUtils:
    """Utility functions for handling data in the Phytovenomics ML platform."""
    
    @staticmethod
    def load_toxin_dataset(file_path: str) -> pd.DataFrame:
        """Load and preprocess the snake toxin dataset.
        
        Args:
            file_path: Path to the toxin dataset CSV file
            
        Returns:
            DataFrame containing the toxin data
        """
        logger.info(f"Loading toxin dataset from {file_path}")
        try:
            if file_path.endswith('.csv'):
                df = pd.read_csv(file_path)
            else:
                raise ValueError(f"Unsupported file format: {file_path}")
                
            logger.info(f"Loaded {len(df)} toxin records")
            return df
        except Exception as e:
            logger.error(f"Error loading toxin dataset: {e}")
            return pd.DataFrame()
    
    @staticmethod
    def load_antibody_dataset(file_path: str) -> Dict[str, str]:
        """Load antibody sequences from a FASTA file.
        
        Args:
            file_path: Path to the antibody FASTA file
            
        Returns:
            Dictionary mapping sequence IDs to sequences
        """
        logger.info(f"Loading antibody dataset from {file_path}")
        sequences = {}
        
        try:
            current_id = None
            current_seq = []
            
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save previous sequence if exists
                        if current_id and current_seq:
                            sequences[current_id] = ''.join(current_seq)
                            
                        # Start new sequence
                        current_id = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
                        
            # Save the last sequence
            if current_id and current_seq:
                sequences[current_id] = ''.join(current_seq)
                
            logger.info(f"Loaded {len(sequences)} antibody sequences")
            return sequences
        except Exception as e:
            logger.error(f"Error loading antibody dataset: {e}")
            return {}
    
    @staticmethod
    def load_binding_data(file_path: str) -> pd.DataFrame:
        """Load antibody-toxin binding affinity data.
        
        Args:
            file_path: Path to the binding data CSV file
            
        Returns:
            DataFrame containing binding data
        """
        logger.info(f"Loading binding data from {file_path}")
        try:
            if file_path.endswith('.csv'):
                df = pd.read_csv(file_path)
            else:
                raise ValueError(f"Unsupported file format: {file_path}")
                
            logger.info(f"Loaded {len(df)} binding data points")
            return df
        except Exception as e:
            logger.error(f"Error loading binding data: {e}")
            return pd.DataFrame()
    
    @staticmethod
    def save_designed_antibodies(antibodies: List[Dict], output_file: str):
        """Save designed antibodies to a JSON file.
        
        Args:
            antibodies: List of antibody dictionaries
            output_file: Path to the output file
        """
        logger.info(f"Saving {len(antibodies)} designed antibodies to {output_file}")
        try:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as f:
                json.dump({"antibodies": antibodies}, f, indent=2)
            logger.info(f"Antibodies saved successfully")
        except Exception as e:
            logger.error(f"Error saving antibodies: {e}")
    
    @staticmethod
    def save_cocktail_design(cocktail: Dict, output_file: str):
        """Save cocktail design to a JSON file.
        
        Args:
            cocktail: Cocktail dictionary
            output_file: Path to the output file
        """
        logger.info(f"Saving cocktail design to {output_file}")
        try:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as f:
                json.dump(cocktail, f, indent=2)
            logger.info(f"Cocktail design saved successfully")
        except Exception as e:
            logger.error(f"Error saving cocktail design: {e}")
    
    @staticmethod
    def generate_summary_report(antibodies: List[Dict], output_file: str):
        """Generate a summary report of designed antibodies.
        
        Args:
            antibodies: List of antibody dictionaries
            output_file: Path to the output report file
        """
        logger.info(f"Generating summary report for {len(antibodies)} antibodies")
        try:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
            # Create summary data
            summary = []
            for ab in antibodies:
                summary_entry = {
                    "id": ab.get("id", "unknown"),
                    "type": ab.get("type", "unknown"),
                    "target": ab.get("target", "unknown"),
                    "binding_energy": ab.get("predicted_binding", {}).get("energy", 0),
                    "affinity_nm": ab.get("predicted_binding", {}).get("affinity_nm", 0),
                    "developability": ab.get("developability", {}).get("overall_score", 0)
                }
                summary.append(summary_entry)
            
            # Write to CSV
            with open(output_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=summary[0].keys())
                writer.writeheader()
                writer.writerows(summary)
                
            logger.info(f"Summary report saved to {output_file}")
        except Exception as e:
            logger.error(f"Error generating summary report: {e}")
            
    @staticmethod
    def convert_to_fasta(antibodies: List[Dict], output_file: str):
        """Convert antibody sequences to FASTA format.
        
        Args:
            antibodies: List of antibody dictionaries
            output_file: Path to the output FASTA file
        """
        logger.info(f"Converting {len(antibodies)} antibodies to FASTA format")
        try:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
            with open(output_file, 'w') as f:
                for ab in antibodies:
                    ab_id = ab.get("id", "unknown")
                    
                    # Write heavy chain
                    if "heavy_chain" in ab and "sequence" in ab["heavy_chain"]:
                        f.write(f">{ab_id}_HC\n{ab['heavy_chain']['sequence']}\n")
                    
                    # Write light chain
                    if "light_chain" in ab and "sequence" in ab["light_chain"]:
                        f.write(f">{ab_id}_LC\n{ab['light_chain']['sequence']}\n")
                        
            logger.info(f"FASTA file saved to {output_file}")
        except Exception as e:
            logger.error(f"Error converting to FASTA: {e}")