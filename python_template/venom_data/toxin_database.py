#!/usr/bin/env python
# venom_data/toxin_database.py

import logging
import pandas as pd
import numpy as np
import os
from typing import Dict, List, Optional, Any, Union
import json
from collections import defaultdict

logger = logging.getLogger("Phytovenomics.ToxinDatabase")

class ToxinDatabase:
    """Database for managing and querying snake toxin data."""
    
    def __init__(self, config: Dict):
        """
        Initialize the ToxinDatabase with configuration.
        
        Args:
            config: Dictionary containing configuration parameters
        """
        self.config = config
        self.toxins = {}
        self.metadata = {}
        self.cache_enabled = config.get("cache_enabled", True)
        self.cache_dir = config.get("cache_dir", ".cache/toxins")
        self.toxin_families = config.get("toxin_families", {})
        self.max_toxins = config.get("max_toxins", 1000)
        
        if self.cache_enabled:
            os.makedirs(self.cache_dir, exist_ok=True)
    
    def load_from_csv(self, file_path: str) -> bool:
        """Load toxin data from a CSV file.
        
        Args:
            file_path: Path to the toxin dataset CSV file
            
        Returns:
            Boolean indicating success
        """
        logger.info(f"Loading toxin data from {file_path}")
        try:
            if not os.path.exists(file_path):
                logger.error(f"Toxin dataset file not found: {file_path}")
                return False
            
            # Read the CSV file
            df = pd.read_csv(file_path)
            
            # Check required columns
            required_columns = ["toxin_id", "sequence", "toxin_family", "source_species"]
            missing_columns = [col for col in required_columns if col not in df.columns]
            
            if missing_columns:
                logger.error(f"Missing required columns: {missing_columns}")
                return False
            
            # Convert to dictionary format
            toxin_count = 0
            for _, row in df.iterrows():
                if toxin_count >= self.max_toxins:
                    logger.warning(f"Reached maximum toxin limit: {self.max_toxins}")
                    break
                    
                toxin_id = row["toxin_id"]
                
                # Create toxin entry
                toxin = {
                    "id": toxin_id,
                    "sequence": row["sequence"],
                    "family": row["toxin_family"],
                    "source_species": row["source_species"]
                }
                
                # Add optional fields if present
                for field in ["molecular_weight", "protein_type", "function", "target", "potency"]:
                    if field in df.columns and not pd.isna(row[field]):
                        toxin[field] = row[field]
                
                self.toxins[toxin_id] = toxin
                toxin_count += 1
            
            # Update metadata
            self.metadata = {
                "source": file_path,
                "toxin_count": len(self.toxins),
                "families": list(set(toxin["family"] for toxin in self.toxins.values()))
            }
            
            logger.info(f"Loaded {len(self.toxins)} toxins from {file_path}")
            
            # Cache data if enabled
            if self.cache_enabled:
                self._cache_toxins()
            
            return True
            
        except Exception as e:
            logger.error(f"Error loading toxin data: {e}")
            return False
    
    def get_toxin_by_id(self, toxin_id: str) -> Optional[Dict]:
        """Get toxin data by ID.
        
        Args:
            toxin_id: Toxin identifier
            
        Returns:
            Dictionary containing toxin data or None if not found
        """
        if toxin_id in self.toxins:
            return self.toxins[toxin_id]
        else:
            logger.warning(f"Toxin ID not found: {toxin_id}")
            return None
    
    def get_toxins_by_family(self, family: str) -> List[Dict]:
        """Get all toxins belonging to a specific family.
        
        Args:
            family: Toxin family name
            
        Returns:
            List of toxin dictionaries
        """
        return [toxin for toxin in self.toxins.values() if toxin["family"] == family]
    
    def get_all_toxins(self) -> Dict[str, Dict]:
        """Get all toxins in the database.
        
        Returns:
            Dictionary mapping toxin IDs to toxin dictionaries
        """
        return self.toxins
    
    def search_toxins(self, query: Dict) -> List[Dict]:
        """Search for toxins matching specific criteria.
        
        Args:
            query: Dictionary containing search criteria
            
        Returns:
            List of matching toxin dictionaries
        """
        results = []
        
        for toxin in self.toxins.values():
            match = True
            
            for key, value in query.items():
                if key not in toxin or toxin[key] != value:
                    match = False
                    break
            
            if match:
                results.append(toxin)
        
        return results
    
    def get_toxin_statistics(self) -> Dict[str, Any]:
        """Get statistics about toxins in the database.
        
        Returns:
            Dictionary containing toxin statistics
        """
        if not self.toxins:
            return {}
        
        # Count toxins by family
        family_counts = defaultdict(int)
        for toxin in self.toxins.values():
            family_counts[toxin["family"]] += 1
        
        # Calculate sequence length statistics
        lengths = [len(toxin["sequence"]) for toxin in self.toxins.values()]
        
        return {
            "total_toxins": len(self.toxins),
            "families": dict(family_counts),
            "sequence_length": {
                "min": min(lengths),
                "max": max(lengths),
                "mean": np.mean(lengths),
                "median": np.median(lengths)
            }
        }
    
    def _cache_toxins(self):
        """Cache toxin data to file."""
        try:
            cache_file = os.path.join(self.cache_dir, "toxin_cache.json")
            cache_data = {
                "metadata": self.metadata,
                "toxins": self.toxins
            }
            
            with open(cache_file, 'w') as f:
                json.dump(cache_data, f)
                
            logger.info(f"Toxin data cached to {cache_file}")
        except Exception as e:
            logger.error(f"Error caching toxin data: {e}")
    
    def load_from_cache(self) -> bool:
        """Load toxin data from cache.
        
        Returns:
            Boolean indicating success
        """
        if not self.cache_enabled:
            logger.warning("Cache is disabled, cannot load from cache")
            return False
            
        try:
            cache_file = os.path.join(self.cache_dir, "toxin_cache.json")
            
            if not os.path.exists(cache_file):
                logger.warning(f"Cache file not found: {cache_file}")
                return False
                
            with open(cache_file, 'r') as f:
                cache_data = json.load(f)
                
            self.metadata = cache_data.get("metadata", {})
            self.toxins = cache_data.get("toxins", {})
            
            logger.info(f"Loaded {len(self.toxins)} toxins from cache")
            return True
            
        except Exception as e:
            logger.error(f"Error loading from cache: {e}")
            return False
    
    def identify_conserved_regions(self, family: str, conservation_threshold: float = 0.7) -> List[Dict]:
        """Identify conserved regions within a toxin family.
        
        Args:
            family: Toxin family name
            conservation_threshold: Minimum conservation score (0-1)
            
        Returns:
            List of dictionaries describing conserved regions
        """
        family_toxins = self.get_toxins_by_family(family)
        
        if len(family_toxins) < 2:
            logger.warning(f"Not enough toxins in family {family} to identify conserved regions")
            return []
            
        # This is a simplified implementation
        # In practice, would use multiple sequence alignment and more sophisticated conservation analysis
        
        # For demonstration, we'll identify regions where at least conservation_threshold fraction of 
        # sequences have the same amino acid at each position
        
        # First, find the length of the shortest sequence
        min_length = min(len(toxin["sequence"]) for toxin in family_toxins)
        
        conserved_regions = []
        current_region = {"start": -1, "end": -1, "sequence": ""}
        
        for i in range(min_length):
            # Count amino acids at this position
            position_aas = [toxin["sequence"][i] for toxin in family_toxins]
            aa_counts = defaultdict(int)
            
            for aa in position_aas:
                aa_counts[aa] += 1
                
            # Find most common amino acid and its frequency
            most_common_aa = max(aa_counts.items(), key=lambda x: x[1])
            conservation = most_common_aa[1] / len(family_toxins)
            
            if conservation >= conservation_threshold:
                # This position is conserved
                if current_region["start"] == -1:
                    # Start a new region
                    current_region = {
                        "start": i,
                        "end": i,
                        "sequence": most_common_aa[0]
                    }
                else:
                    # Extend current region
                    current_region["end"] = i
                    current_region["sequence"] += most_common_aa[0]
            elif current_region["start"] != -1:
                # End of conserved region
                if len(current_region["sequence"]) >= 5:  # Minimum length for a region
                    conserved_regions.append({
                        "family": family,
                        "start": current_region["start"],
                        "end": current_region["end"],
                        "length": current_region["end"] - current_region["start"] + 1,
                        "sequence": current_region["sequence"],
                        "conservation": conservation_threshold
                    })
                
                # Reset current region
                current_region = {"start": -1, "end": -1, "sequence": ""}
        
        # Check if we have a final region to add
        if current_region["start"] != -1 and len(current_region["sequence"]) >= 5:
            conserved_regions.append({
                "family": family,
                "start": current_region["start"],
                "end": current_region["end"],
                "length": current_region["end"] - current_region["start"] + 1,
                "sequence": current_region["sequence"],
                "conservation": conservation_threshold
            })
        
        return conserved_regions