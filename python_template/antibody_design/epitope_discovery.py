#!/usr/bin/env python
# antibody_design/epitope_discovery.py

import logging
import random
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Dict, List, Optional, Tuple, Set
from utils.model_utils import EpitopePredictor, ModelUtils
from pathlib import Path
import pandas as pd
import joblib

logger = logging.getLogger("Phytovenomics.EpitopeDiscovery")

class EpitopeDiscovery:
    """Discovers potential epitopes on toxin proteins for antibody targeting."""
    
    def __init__(self, config: Dict):
        """
        Initialize the EpitopeDiscovery with configuration.
        
        Args:
            config: Dictionary containing configuration parameters
        """
        self.config = config
        self.min_epitope_length = config.get("min_epitope_length", 8)
        self.max_epitope_length = config.get("max_epitope_length", 20)
        self.ideal_epitope_length = config.get("ideal_epitope_length", 15)
        
        # ML model configuration
        self.use_ml_model = config.get("use_ml_model", True)
        self.model_path = config.get("epitope_model_path", "models/best_epitope_predictor.pt")
        
        # Device for ML models
        self.device = ModelUtils.get_device() if torch.cuda.is_available() else torch.device("cpu")
        
        # Initialize ML model if enabled
        if self.use_ml_model:
            self._initialize_ml_model()
        
        # Weights for scoring epitope candidates
        self.scoring_weights = config.get("scoring_weights", {
            "ml_prediction": 0.4,  # New ML-based prediction weight
            "accessibility": 0.2,
            "hydrophilicity": 0.15,
            "conservation": 0.15,
            "immunogenicity": 0.1
        })
        
        # Thresholds
        self.conservation_threshold = config.get("conservation_threshold", 0.6)
        self.accessibility_threshold = config.get("accessibility_threshold", 0.5)
        self.epitope_overlap_threshold = config.get("epitope_overlap_threshold", 5)
        
        # Amino acid properties
        self.hydrophilicity_scale = {
            "A": -0.5, "R": 3.0, "N": 0.2, "D": 3.0, "C": -1.0,
            "Q": 0.2, "E": 3.0, "G": 0.0, "H": -0.5, "I": -1.8,
            "L": -1.8, "K": 3.0, "M": -1.3, "F": -2.5, "P": 0.0,
            "S": 0.3, "T": -0.4, "W": -3.4, "Y": -2.3, "V": -1.5
        }
        
        # Simple immunogenicity scale (higher values = more immunogenic)
        self.immunogenicity_scale = {
            "A": 0.3, "R": 0.7, "N": 0.5, "D": 0.7, "C": 0.5,
            "Q": 0.5, "E": 0.7, "G": 0.3, "H": 0.8, "I": 0.5,
            "L": 0.5, "K": 0.7, "M": 0.5, "F": 0.8, "P": 0.3,
            "S": 0.4, "T": 0.4, "W": 0.9, "Y": 0.9, "V": 0.5
        }
        
        # Amino acid encoding for ML model
        self.aa_to_idx = {
            'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 
            'K': 9, 'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 
            'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20, 'X': 0, '-': 0, '*': 0
        }
        
        # Feature extraction parameters
        self.window_size = config.get("epitope_window_size", 21)  # Window size for sequence context
        
    def _initialize_ml_model(self):
        """Initialize the ML model for epitope prediction."""
        try:
            # Create model directory if it doesn't exist
            model_dir = Path(self.model_path).parent
            model_dir.mkdir(exist_ok=True, parents=True)
            
            # Initialize model
            self.model = EpitopePredictor(embed_dim=512, hidden_dim=256, num_layers=2)
            
            # Check if model file exists, if not create a basic model
            model_file = Path(self.model_path)
            if model_file.exists():
                logger.info(f"Loading epitope prediction model from {self.model_path}")
                state_dict = torch.load(self.model_path, map_location=self.device)
                self.model.load_state_dict(state_dict)
            else:
                logger.warning(f"Model file {self.model_path} not found. Using untrained model.")
                # Save the untrained model for future reference
                torch.save(self.model.state_dict(), self.model_path)
                
            self.model.to(self.device)
            self.model.eval()  # Set model to evaluation mode
            logger.info("ML model for epitope prediction initialized")
            
            # Initialize feature scaler for sequence encoding
            scaler_path = Path(self.model_path).parent / "epitope_feature_scaler.joblib"
            if scaler_path.exists():
                self.feature_scaler = joblib.load(scaler_path)
                logger.info(f"Loaded feature scaler from {scaler_path}")
            else:
                self.feature_scaler = None
                logger.warning(f"Feature scaler not found at {scaler_path}. Will normalize on-the-fly.")
                
        except Exception as e:
            logger.error(f"Error initializing ML model: {e}")
            logger.info("Falling back to rule-based epitope prediction")
            self.use_ml_model = False
    
    def identify_epitopes(self, toxin: Dict) -> List[Dict]:
        """
        Identify potential epitopes on a toxin protein using ML model when available.
        
        Args:
            toxin: Dictionary containing toxin data
            
        Returns:
            List of dictionaries describing identified epitopes
        """
        logger.info(f"Identifying epitopes for toxin {toxin.get('id')}")
        
        sequence = toxin.get("sequence", "")
        if not sequence:
            logger.warning(f"Toxin {toxin.get('id')} has no sequence data")
            return []
        
        # Generate epitope candidates
        candidates = self._generate_epitope_candidates(sequence)
        
        # Score candidates using both rule-based and ML-based approaches
        if self.use_ml_model and hasattr(self, 'model'):
            logger.info(f"Using ML-enhanced epitope prediction for toxin {toxin.get('id')}")
            try:
                # Enhance with ML predictions
                scored_candidates = self._score_epitope_candidates(candidates, toxin)
                logger.info(f"ML model successfully scored {len(scored_candidates)} epitope candidates")
            except Exception as e:
                logger.error(f"Error in ML epitope prediction: {e}")
                logger.info("Falling back to rule-based prediction")
                # Reset ML flag for this run
                orig_ml_flag = self.use_ml_model
                self.use_ml_model = False
                scored_candidates = self._score_epitope_candidates(candidates, toxin)
                self.use_ml_model = orig_ml_flag
        else:
            # Use only rule-based approach
            logger.info(f"Using rule-based epitope prediction for toxin {toxin.get('id')}")
            scored_candidates = self._score_epitope_candidates(candidates, toxin)
        
        # Sort by score (descending)
        scored_candidates.sort(key=lambda x: x["score"], reverse=True)
        
        # Filter overlapping epitopes
        filtered_epitopes = self._filter_overlapping_epitopes(scored_candidates)
        
        # Add metadata to epitopes
        for i, epitope in enumerate(filtered_epitopes):
            epitope_id = f"{toxin.get('id')}_ep{i+1}"
            epitope.update({
                "id": epitope_id,
                "toxin_id": toxin.get("id"),
                "toxin_family": toxin.get("family"),
                "source_species": toxin.get("source_species"),
                "prediction_method": "ML-enhanced" if self.use_ml_model else "Rule-based"
            })
        
        logger.info(f"Identified {len(filtered_epitopes)} epitopes for toxin {toxin.get('id')}")
        return filtered_epitopes
        

            
    def predict_epitopes_with_ml(self, toxin_sequence: str) -> List[Dict]:
        """
        Perform a complete end-to-end ML-based epitope prediction on a toxin sequence.
        
        Args:
            toxin_sequence: Full amino acid sequence of the toxin
            
        Returns:
            List of predicted epitope regions with scores
        """
        if not self.use_ml_model or not hasattr(self, 'model'):
            logger.warning("ML model not enabled or initialized. Using rule-based prediction.")
            return self.identify_epitopes({"id": "unknown", "sequence": toxin_sequence})
            
        logger.info(f"Performing end-to-end ML epitope prediction on sequence of length {len(toxin_sequence)}")
        
        try:
            # Ensure model is in evaluation mode
            self.model.eval()
            
            # Process the full sequence with sliding window
            windows = []
            positions = []
            window_size = self.window_size
            half_window = window_size // 2
            
            # Pad sequence for edge cases
            padded_seq = "X" * half_window + toxin_sequence + "X" * half_window
            
            # Generate windows for each position
            for i in range(len(toxin_sequence)):
                window = padded_seq[i:i + window_size]
                if len(window) == window_size:  # Ensure full window
                    windows.append(window)
                    positions.append(i)
            
            # Convert windows to tensors
            window_tensors = []
            batch_size = 32  # Process in batches
            
            for i in range(0, len(windows), batch_size):
                batch_windows = windows[i:i+batch_size]
                batch_embeddings = []
                
                # Process all windows in the batch at once with consistent dimensions
                for window in batch_windows:
                    # Create simple embedding for window - ensure all windows have EXACTLY same size
                    embedding = torch.zeros((1, window_size, 512))
                    
                    # Only process up to window_size characters
                    for j, aa in enumerate(window[:window_size]):
                        # Fill in with amino acid properties
                        hydro = self.hydrophilicity_scale.get(aa, 0)
                        immuno = self.immunogenicity_scale.get(aa, 0.5)
                        aa_idx = self.aa_to_idx.get(aa.upper(), 0)
                        
                        embedding[0, j, 0] = hydro
                        embedding[0, j, 1] = immuno
                        embedding[0, j, aa_idx % 512] = 1.0
                        
                    batch_embeddings.append(embedding)
                
                if batch_embeddings:
                    window_tensors.extend(batch_embeddings)
            
            # Make predictions
            predictions = []
            with torch.no_grad():
                for tensor in window_tensors:
                    tensor = tensor.to(self.device)
                    pred = self.model(tensor)
                    predictions.append(pred.mean().item())  # Mean across window
            
            # Combine predictions with positions
            pos_scores = [(pos, score) for pos, score in zip(positions, predictions)]
            
            # Find continuous regions with high scores
            epitope_regions = []
            current_region = []
            threshold = 0.6  # Threshold for positive prediction
            
            for pos, score in pos_scores:
                if score >= threshold:
                    if not current_region:
                        current_region = [(pos, score)]
                    elif pos == current_region[-1][0] + 1:  # Adjacent position
                        current_region.append((pos, score))
                    else:  # Gap in positions
                        # Save completed region if long enough
                        if len(current_region) >= self.min_epitope_length:
                            epitope_regions.append(current_region)
                        # Start new region
                        current_region = [(pos, score)]
                else:
                    # End of a region
                    if current_region and len(current_region) >= self.min_epitope_length:
                        epitope_regions.append(current_region)
                        current_region = []
                    else:
                        current_region = []
            
            # Check last region
            if current_region and len(current_region) >= self.min_epitope_length:
                epitope_regions.append(current_region)
            
            # Convert regions to epitope dictionaries
            epitopes = []
            for i, region in enumerate(epitope_regions):
                start_pos = region[0][0]
                end_pos = region[-1][0]
                avg_score = sum(s for _, s in region) / len(region)
                sequence = toxin_sequence[start_pos:end_pos+1]
                
                epitopes.append({
                    "id": f"ml_epitope_{i+1}",
                    "start": start_pos,
                    "end": end_pos,
                    "length": end_pos - start_pos + 1,
                    "sequence": sequence,
                    "ml_score": round(avg_score, 2),
                    "prediction_method": "ML-direct"
                })
            
            logger.info(f"ML-direct prediction identified {len(epitopes)} epitope regions")
            return epitopes
            
        except Exception as e:
            logger.error(f"Error in end-to-end ML epitope prediction: {e}")
            logger.info("Falling back to standard epitope identification")
            return self.identify_epitopes({"id": "unknown", "sequence": toxin_sequence})
    
    def _generate_epitope_candidates(self, sequence: str) -> List[Dict]:
        """Generate potential epitope candidates from a protein sequence."""
        candidates = []
        
        # Generate all possible epitopes within length constraints
        seq_len = len(sequence)
        for length in range(self.min_epitope_length, min(self.max_epitope_length + 1, seq_len + 1)):
            for start in range(seq_len - length + 1):
                end = start + length - 1
                epitope_seq = sequence[start:start+length]
                
                candidates.append({
                    "start": start,
                    "end": end,
                    "length": length,
                    "sequence": epitope_seq
                })
        
        return candidates
    
    def _score_epitope_candidates(self, candidates: List[Dict], toxin: Dict) -> List[Dict]:
        """Score epitope candidates based on multiple criteria including ML predictions."""
        scored_candidates = []
        
        # Prepare all candidates for batch ML prediction if ML model is enabled
        ml_predictions = {}
        if self.use_ml_model:
            try:
                ml_predictions = self._predict_epitopes_with_ml(toxin["sequence"], candidates)
                logger.info(f"Generated ML predictions for {len(ml_predictions)} epitope candidates")
            except Exception as e:
                logger.error(f"Error using ML model for epitope prediction: {e}")
                logger.info("Falling back to rule-based scoring only")
        
        for candidate in candidates:
            # Calculate various scores
            accessibility = self._calculate_accessibility_score(candidate["sequence"])
            hydrophilicity = self._calculate_hydrophilicity_score(candidate["sequence"])
            conservation = self._calculate_conservation_score(candidate, toxin)
            immunogenicity = self._calculate_immunogenicity_score(candidate["sequence"])
            
            # Get ML prediction score if available
            ml_score = ml_predictions.get(candidate["sequence"], 0.5)  # Default to neutral if no ML prediction
            
            # Apply length penalty (prefer epitopes closer to ideal length)
            length_factor = 1.0 - abs(candidate["length"] - self.ideal_epitope_length) / 10.0
            length_factor = max(0.7, min(1.0, length_factor))  # Limit penalty
            
            # Calculate weighted score
            weights = self.scoring_weights
            score = (
                (weights.get("ml_prediction", 0.0) * ml_score if self.use_ml_model else 0.0) +
                weights.get("accessibility", 0.3) * accessibility +
                weights.get("hydrophilicity", 0.2) * hydrophilicity +
                weights.get("conservation", 0.25) * conservation +
                weights.get("immunogenicity", 0.25) * immunogenicity
            ) * length_factor
            
            # Add scores to candidate
            candidate_with_scores = candidate.copy()
            candidate_with_scores.update({
                "accessibility": round(accessibility, 2),
                "hydrophilicity": round(hydrophilicity, 2),
                "conservation": round(conservation, 2),
                "predicted_immunogenicity": round(immunogenicity, 2),
                "ml_score": round(ml_score, 2) if self.use_ml_model else None,
                "score": round(score, 2)
            })
            
            # Apply basic filtering
            if accessibility >= self.accessibility_threshold and conservation >= self.conservation_threshold:
                scored_candidates.append(candidate_with_scores)
        
        return scored_candidates
    
    def _predict_epitopes_with_ml(self, toxin_sequence: str, candidates: List[Dict]) -> Dict[str, float]:
        """Use the ML model to predict epitope scores for candidate sequences.
        
        Args:
            toxin_sequence: Full toxin sequence
            candidates: List of epitope candidates
            
        Returns:
            Dictionary mapping epitope sequences to their ML prediction scores (0-1)
        """
        # Ensure model is in evaluation mode
        self.model.eval()
        
        # We'll use a simpler approach to avoid tensor size mismatches
        # Instead of using the encoded full sequence, we'll encode each epitope candidate directly
        predictions = {}
        batch_size = 16  # Process in batches to avoid memory issues
        
        try:
            with torch.no_grad():
                for i in range(0, len(candidates), batch_size):
                    batch = candidates[i:i+batch_size]
                    batch_scores = []
                    
                    for candidate in batch:
                        # Get epitope sequence
                        epitope_seq = candidate["sequence"]
                        
                        # Create a fixed-size embedding for each epitope
                        # Pad or truncate to 8 amino acids (the size expected by the model)
                        fixed_length = 8  # Use a consistent length for all inputs
                        
                        if len(epitope_seq) >= fixed_length:
                            # Take first fixed_length amino acids
                            proc_seq = epitope_seq[:fixed_length]
                        else:
                            # Pad with X to reach fixed_length
                            proc_seq = epitope_seq + "X" * (fixed_length - len(epitope_seq))
                        
                        # Create simple embedding
                        embedding = torch.zeros((1, fixed_length, 512))
                        
                        # Fill embedding with amino acid properties
                        for j, aa in enumerate(proc_seq):
                            hydro = self.hydrophilicity_scale.get(aa, 0)
                            immuno = self.immunogenicity_scale.get(aa, 0.5)
                            aa_idx = self.aa_to_idx.get(aa.upper(), 0)
                            
                            embedding[0, j, 0] = hydro
                            embedding[0, j, 1] = immuno
                            embedding[0, j, aa_idx % 512] = 1.0
                        
                        # Get prediction score for this epitope
                        embedding = embedding.to(self.device)
                        pred = self.model(embedding)
                        score = pred.mean().item()
                        
                        # Ensure score is in 0-1 range
                        score = max(0.0, min(1.0, score))
                        batch_scores.append(score)
                        predictions[epitope_seq] = score
            
            return predictions
            
        except Exception as e:
            logger.error(f"Error in ML prediction: {e}")
            # Return empty predictions dict on error
            return {}
        
        return predictions
    
    def _encode_sequence(self, sequence: str) -> torch.Tensor:
        """Convert a protein sequence to a tensor of embeddings for ML model input.
        
        Args:
            sequence: Amino acid sequence
            
        Returns:
            Tensor of shape [1, seq_len, embed_dim]
        """
        # Convert sequence to indices
        indices = [self.aa_to_idx.get(aa.upper(), 0) for aa in sequence]
        
        # Create tensor of correct shape for model [1, seq_len]
        seq_tensor = torch.tensor([indices], dtype=torch.long)
        
        # Create a simple embedding based on amino acid properties
        # This is a placeholder for actual learned embeddings in a trained model
        embed_dim = 512
        embedding = torch.zeros((1, len(sequence), embed_dim))
        
        # Create simple embeddings based on amino acid properties
        for i, aa in enumerate(sequence):
            # Hydrophilicity contribution
            hydro = self.hydrophilicity_scale.get(aa, 0)
            # Immunogenicity contribution
            immuno = self.immunogenicity_scale.get(aa, 0.5)
            # One-hot encoding component
            aa_idx = self.aa_to_idx.get(aa.upper(), 0)
            
            # Create a sparse embedding
            embedding[0, i, 0] = hydro  # First dimension: hydrophilicity
            embedding[0, i, 1] = immuno  # Second dimension: immunogenicity
            embedding[0, i, aa_idx % embed_dim] = 1.0  # Sparse one-hot component
        
        return embedding
    
    def _filter_overlapping_epitopes(self, epitopes: List[Dict]) -> List[Dict]:
        """Filter out overlapping epitopes, keeping higher-scoring ones."""
        if not epitopes:
            return []
            
        # Sort by score (descending)
        sorted_epitopes = sorted(epitopes, key=lambda x: x["score"], reverse=True)
        
        filtered_epitopes = [sorted_epitopes[0]]
        
        for epitope in sorted_epitopes[1:]:
            # Check if this epitope overlaps significantly with any already selected
            overlapping = False
            for selected in filtered_epitopes:
                if self._check_overlap(epitope, selected):
                    overlapping = True
                    break
            
            if not overlapping:
                filtered_epitopes.append(epitope)
        
        return filtered_epitopes
    
    def _check_overlap(self, epitope1: Dict, epitope2: Dict) -> bool:
        """Check if two epitopes overlap by more than the threshold."""
        # Calculate overlap
        start1, end1 = epitope1["start"], epitope1["end"]
        start2, end2 = epitope2["start"], epitope2["end"]
        
        # Check if there is any overlap
        if end1 < start2 or end2 < start1:
            return False
        
        # Calculate overlap length
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        overlap_length = overlap_end - overlap_start + 1
        
        return overlap_length >= self.epitope_overlap_threshold
    
    def _calculate_accessibility_score(self, sequence: str) -> float:
        """Calculate surface accessibility score for an epitope sequence."""
        # This is a simplified model
        # In a real implementation, would use structural prediction or experimental data
        
        # Amino acids likely to be on the surface (hydrophilic, charged)
        surface_prone = "DEKHRGNQSTP"
        
        # Count surface-prone residues
        surface_count = sum(1 for aa in sequence if aa in surface_prone)
        
        # Calculate score (0-1)
        return surface_count / len(sequence)
    
    def _calculate_hydrophilicity_score(self, sequence: str) -> float:
        """Calculate hydrophilicity score for an epitope sequence."""
        # Calculate average hydrophilicity
        total_hydrophilicity = sum(self.hydrophilicity_scale.get(aa, 0) for aa in sequence)
        avg_hydrophilicity = total_hydrophilicity / len(sequence)
        
        # Normalize to 0-1 scale (from original -3.4 to 3.0 scale)
        normalized_score = (avg_hydrophilicity + 3.4) / 6.4
        
        return min(1.0, max(0.0, normalized_score))  # Ensure 0-1 range
    
    def _calculate_conservation_score(self, candidate: Dict, toxin: Dict) -> float:
        """Calculate conservation score for an epitope candidate."""
        # In a real implementation, would compare with aligned sequences from the same family
        # For this simplified version, we'll use a synthetic conservation score
        
        # For demonstration, assume some regions are more conserved
        # In reality, would be derived from multiple sequence alignment
        sequence = candidate["sequence"]
        start = candidate["start"]
        
        # Simplified model: certain positions in the toxin are more conserved
        # These would normally come from analysis of multiple toxins in the same family
        toxin_id = toxin.get("id", "")
        toxin_family = toxin.get("family", "")
        
        # Seed random generator with toxin ID for reproducibility but varying results
        # This simulates different conservation patterns for different toxins
        seed_value = sum(ord(c) for c in toxin_id) if toxin_id else 42
        random.seed(seed_value)
        
        # Generate synthetic conservation pattern
        # In a real system, this would be based on sequence alignment data
        sequence_length = len(toxin.get("sequence", ""))
        synthetic_conservation = [random.uniform(0.3, 1.0) for _ in range(sequence_length)]
        
        # Make some families more conserved than others
        family_modifier = 1.0
        if "three_finger" in toxin_family.lower():
            family_modifier = 1.2  # More conserved
        elif "phospholipase" in toxin_family.lower():
            family_modifier = 0.9  # Less conserved
        
        # Calculate average conservation for this epitope
        epitope_conservation = sum(synthetic_conservation[start+i] for i in range(len(sequence))) / len(sequence)
        epitope_conservation = min(1.0, epitope_conservation * family_modifier)
        
        return epitope_conservation
    
    def _calculate_immunogenicity_score(self, sequence: str) -> float:
        """Calculate predicted immunogenicity score for an epitope sequence."""
        # Calculate average immunogenicity based on amino acid composition
        total_immunogenicity = sum(self.immunogenicity_scale.get(aa, 0.5) for aa in sequence)
        avg_immunogenicity = total_immunogenicity / len(sequence)
        
        # Look for common immunogenic patterns
        patterns = [
            "YY", "FF", "WW", "KR", "RK", "DE", "ED",  # Adjacent strong residues
            "YFW", "KRH", "DE"  # Clusters of similar properties
        ]
        
        pattern_bonus = 0.0
        for pattern in patterns:
            if pattern in sequence:
                pattern_bonus += 0.05  # Small bonus for each immunogenic pattern
        
        # Combine base score with pattern bonus
        immunogenicity = avg_immunogenicity + pattern_bonus
        
        # Ensure 0-1 range
        return min(1.0, max(0.0, immunogenicity))
    
    def rank_epitopes_for_antibody_design(self, epitopes: List[Dict]) -> List[Dict]:
        """Rank epitopes by their suitability for antibody design."""
        # Create a copy to avoid modifying original
        ranked_epitopes = [epitope.copy() for epitope in epitopes]
        
        # Calculate design suitability score
        for epitope in ranked_epitopes:
            # Factors that make an epitope good for antibody design:
            # - High immunogenicity
            # - Good accessibility
            # - Appropriate length
            # - Conservation (for broadly neutralizing antibodies)
            
            immunogenicity = epitope.get("predicted_immunogenicity", 0)
            accessibility = epitope.get("accessibility", 0)
            conservation = epitope.get("conservation", 0)
            
            # Length factor (prefer epitopes closer to ideal length)
            length = epitope.get("length", 0)
            length_factor = 1.0 - abs(length - self.ideal_epitope_length) / 10.0
            length_factor = max(0.7, min(1.0, length_factor))  # Limit penalty
            
            # Calculate design score
            design_score = (
                0.4 * immunogenicity +
                0.3 * accessibility +
                0.2 * conservation +
                0.1 * length_factor
            )
            
            epitope["design_score"] = round(design_score, 2)
        
        # Sort by design score (descending)
        ranked_epitopes.sort(key=lambda x: x.get("design_score", 0), reverse=True)
        
        return ranked_epitopes
    
    def identify_conserved_epitopes_across_family(self, epitopes_by_toxin: Dict[str, List[Dict]]) -> List[Dict]:
        """Identify epitopes that are conserved across a toxin family."""
        # Group epitopes by family
        family_epitopes = {}
        
        for toxin_id, epitopes in epitopes_by_toxin.items():
            for epitope in epitopes:
                family = epitope.get("toxin_family")
                if family:
                    if family not in family_epitopes:
                        family_epitopes[family] = []
                    family_epitopes[family].append(epitope)
        
        # Find conserved epitopes within each family
        conserved_epitopes = []
        
        for family, epitopes in family_epitopes.items():
            # Need multiple toxins to identify conservation
            if len(set(ep["toxin_id"] for ep in epitopes)) < 2:
                continue
                
            # Group epitopes by sequence similarity
            similar_groups = self._group_similar_epitopes(epitopes)
            
            # For each group with epitopes from multiple toxins, create a conserved epitope
            for group in similar_groups:
                toxin_ids = set(ep["toxin_id"] for ep in group)
                
                if len(toxin_ids) >= 2:  # Epitope appears in at least 2 toxins
                    # Create a representative epitope
                    rep_epitope = group[0].copy()
                    
                    # Update with conservation information
                    rep_epitope.update({
                        "id": f"{family}_conserved_{len(conserved_epitopes)+1}",
                        "conserved_across": list(toxin_ids),
                        "conservation_level": len(toxin_ids) / len(set(ep["toxin_id"] for ep in epitopes)),
                        "family_coverage": len(toxin_ids) / len(set(ep["toxin_id"] for ep in family_epitopes[family])),
                        "is_conserved": True
                    })
                    
                    conserved_epitopes.append(rep_epitope)
        
        # Sort by conservation level (descending)
        conserved_epitopes.sort(key=lambda x: x.get("conservation_level", 0), reverse=True)
        
        return conserved_epitopes
    
    def _group_similar_epitopes(self, epitopes: List[Dict]) -> List[List[Dict]]:
        """Group epitopes by sequence similarity."""
        groups = []
        assigned = set()
        
        for i, ep1 in enumerate(epitopes):
            if i in assigned:
                continue
                
            # Start a new group
            group = [ep1]
            assigned.add(i)
            
            # Find similar epitopes
            for j, ep2 in enumerate(epitopes):
                if j != i and j not in assigned:
                    if self._sequence_similarity(ep1["sequence"], ep2["sequence"]) >= 0.7:  # 70% similarity threshold
                        group.append(ep2)
                        assigned.add(j)
            
            groups.append(group)
        
        return groups
    
    def _sequence_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence similarity between two epitope sequences."""
        # Use simple sequence alignment approach
        # For more accurate results, would use proper alignment algorithm
        
        # If lengths are too different, similarity is lower
        len_diff = abs(len(seq1) - len(seq2))
        if len_diff / max(len(seq1), len(seq2)) > 0.3:  # More than 30% length difference
            return 0.0
        
        # Find subsequence matches
        matches = 0
        shorter = seq1 if len(seq1) <= len(seq2) else seq2
        longer = seq2 if len(seq1) <= len(seq2) else seq1
        
        # Try different alignments and take the best
        best_match = 0
        
        # Slide shorter sequence across longer one
        for start in range(len(longer) - len(shorter) + 1):
            current_matches = sum(a == b for a, b in zip(shorter, longer[start:start+len(shorter)]))
            best_match = max(best_match, current_matches)
        
        return best_match / len(shorter)