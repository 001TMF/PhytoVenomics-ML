#!/usr/bin/env python
# antibody_design/esm2_cdr_masking.py

import os
import re
import time
import yaml
import torch
import logging
import random
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Union, Any
from pathlib import Path
from transformers import AutoTokenizer, AutoModel, AutoModelForMaskedLM
from torch.nn import functional as F

from antibody_design.cdr_processor import CDRProcessor
from utils.training_utils import setup_device, set_seed

logger = logging.getLogger("Phytovenomics.ESM2CDR")

@dataclass
class ESM2Config:
    """
    Configuration for ESM-2 CDR masking model
    """
    model: Dict[str, Any] = field(default_factory=lambda: {
        "esm2_model": "esm2_t33_650M_UR50D",
        "esm2_pretrained": True,
        "use_esm_embeddings": True,
        "freeze_pretrained": False,
        "embedding_dim": 1280,
        "hidden_dim": 512,
        "num_attention_heads": 20,
        "num_hidden_layers": 2,
        "dropout_rate": 0.1,
        "max_sequence_length": 512
    })
    
    training: Dict[str, Any] = field(default_factory=lambda: {
        "learning_rate": 5.0e-5,
        "warmup_steps": 100,
        "weight_decay": 0.01,
        "batch_size": 8,
        "epochs": 20,
        "early_stopping_patience": 3,
        "gradient_accumulation_steps": 4,
        "max_grad_norm": 1.0,
        "use_fp16": True,
        "seed": 42
    })
    
    cdr_masking: Dict[str, Any] = field(default_factory=lambda: {
        "mask_prob": 0.8,
        "cdr_regions": ["CDRH1", "CDRH2", "CDRH3", "CDRL1", "CDRL2", "CDRL3"],
        "cdr_weighting": True,
        "cdr_weight_factor": 2.0,
        "framework_mask_prob": 0.15,
        "whole_cdr_masking": True
    })
    
    antibody_generation: Dict[str, Any] = field(default_factory=lambda: {
        "num_beams": 5,
        "top_k": 50,
        "top_p": 0.95,
        "temperature": 0.7,
        "repetition_penalty": 1.2,
        "no_repeat_ngram_size": 3,
        "max_generation_attempts": 10,
        "min_cdr_length": 3,
        "max_cdr_length": 20
    })
    
    binding_prediction: Dict[str, Any] = field(default_factory=lambda: {
        "use_cross_attention": True,
        "toxin_embedding_dim": 1280,
        "binding_threshold": 0.5,
        "epitope_attention": True,
        "use_binding_data": True
    })
    
    dataset: Dict[str, Any] = field(default_factory=lambda: {
        "train_split": 0.8,
        "val_split": 0.1,
        "test_split": 0.1,
        "augmentation_factor": 2,
        "use_binding_pairs": True,
        "max_training_examples": 10000
    })
    
    checkpoint: Dict[str, Any] = field(default_factory=lambda: {
        "save_dir": "checkpoints/esm2_cdr",
        "save_steps": 1000,
        "save_total_limit": 5,
        "logging_steps": 100,
        "eval_steps": 500,
        "use_wandb": False
    })
    
    device: str = "auto"
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> "ESM2Config":
        """Create config from dictionary"""
        config = cls()
        for key, value in config_dict.items():
            if hasattr(config, key):
                setattr(config, key, value)
        return config
    
    @classmethod
    def from_yaml(cls, config_path: str) -> "ESM2Config":
        """Load config from YAML file"""
        with open(config_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        return cls.from_dict(config_dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary"""
        return {
            "model": self.model,
            "training": self.training,
            "cdr_masking": self.cdr_masking,
            "antibody_generation": self.antibody_generation,
            "binding_prediction": self.binding_prediction,
            "dataset": self.dataset,
            "checkpoint": self.checkpoint,
            "device": self.device
        }
    
    def save_to_yaml(self, config_path: str) -> None:
        """Save config to YAML file"""
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        with open(config_path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False)


class ESM2CDRModel:
    """
    ESM-2 model for CDR prediction and antibody design.
    
    This class implements:
    1. CDR masking and prediction
    2. Antibody sequence generation
    3. Toxin-antibody binding prediction
    """
    
    def __init__(self, config: ESM2Config = None):
        """
        Initialize ESM-2 CDR Model.
        
        Args:
            config: ESM-2 configuration
        """
        # Initialize configuration
        self.config = config if config is not None else ESM2Config()
        
        # Set random seed
        set_seed(self.config.training["seed"])
        
        # Set up device
        self.device = setup_device(self.config.device)
        logger.info(f"Using device: {self.device}")
        
        # Initialize CDR processor
        self.cdr_processor = CDRProcessor()
        
        # Load ESM-2 model and tokenizer
        logger.info(f"Loading ESM-2 model: {self.config.model['esm2_model']}")
        try:
            self.tokenizer = AutoTokenizer.from_pretrained(self.config.model["esm2_model"])
            
            # For prediction tasks, we need the masked language model
            self.model = AutoModelForMaskedLM.from_pretrained(
                self.config.model["esm2_model"], 
                local_files_only=False
            )
            
            # For embedding tasks, we use the base model
            self.base_model = AutoModel.from_pretrained(
                self.config.model["esm2_model"],
                local_files_only=False
            )
            
            if not self.config.model["freeze_pretrained"]:
                logger.info("Fine-tuning mode: ESM-2 model weights will be updated during training")
            else:
                logger.info("Feature extraction mode: ESM-2 model weights are frozen")
                for param in self.model.parameters():
                    param.requires_grad = False
                for param in self.base_model.parameters():
                    param.requires_grad = False
            
            # Move models to device
            self.model.to(self.device)
            self.base_model.to(self.device)
            
            logger.info("ESM-2 model and tokenizer loaded successfully")
        except Exception as e:
            logger.error(f"Error loading ESM-2 model: {e}")
            raise RuntimeError(f"Failed to load ESM-2 model: {e}")
        
        # Initialize head for binding prediction if needed
        if hasattr(self.config, "binding_prediction") and self.config.binding_prediction["use_cross_attention"]:
            self._init_binding_prediction_head()
    
    def _init_binding_prediction_head(self):
        """
        Initialize the binding prediction head with cross-attention.
        """
        # Create simple MLP for binding prediction
        try:
            import torch.nn as nn
            
            emb_dim = self.config.model["embedding_dim"]
            hidden_dim = self.config.model["hidden_dim"]
            
            # Cross-attention modules
            self.query_proj = nn.Linear(emb_dim, hidden_dim).to(self.device)
            self.key_proj = nn.Linear(emb_dim, hidden_dim).to(self.device)
            self.value_proj = nn.Linear(emb_dim, hidden_dim).to(self.device)
            
            # Binding prediction layers
            self.binding_layers = nn.Sequential(
                nn.Linear(hidden_dim * 2, hidden_dim),
                nn.ReLU(),
                nn.Dropout(self.config.model["dropout_rate"]),
                nn.Linear(hidden_dim, 1),
                nn.Sigmoid()
            ).to(self.device)
            
            logger.info("Binding prediction head initialized")
        except Exception as e:
            logger.error(f"Error initializing binding prediction head: {e}")
            raise RuntimeError(f"Failed to initialize binding prediction head: {e}")
    
    def encode_sequence(self, sequence: str) -> torch.Tensor:
        """
        Encode a protein sequence using ESM-2 model.
        
        Args:
            sequence: Amino acid sequence
            
        Returns:
            Tensor of sequence embeddings
        """
        # Clean sequence
        sequence = re.sub(r"[^A-Za-z]", "", sequence).upper()
        
        # Tokenize sequence
        inputs = self.tokenizer(sequence, return_tensors="pt")
        inputs = {k: v.to(self.device) for k, v in inputs.items()}
        
        # Get embeddings from base model
        with torch.no_grad():
            outputs = self.base_model(**inputs)
            # Get mean of token embeddings (excluding special tokens)
            embeddings = outputs.last_hidden_state[:, 1:-1, :].mean(dim=1)
            
        return embeddings
    
    def mask_cdr(self, sequence: str, cdr_name: str) -> Optional[str]:
        """
        Mask a specific CDR in the antibody sequence.
        
        Args:
            sequence: Antibody amino acid sequence
            cdr_name: Name of the CDR to mask (e.g., "CDRH3")
            
        Returns:
            Sequence with the specified CDR masked, or None if CDR not found
        """
        return self.cdr_processor.mask_cdr_in_sequence(sequence, cdr_name)
    
    def mask_cdrs(self, sequence: str, cdr_names: List[str] = None) -> Optional[str]:
        """
        Mask multiple CDRs in the antibody sequence.
        
        Args:
            sequence: Antibody amino acid sequence
            cdr_names: List of CDR names to mask (e.g., ["CDRH3", "CDRL3"])
            
        Returns:
            Sequence with specified CDRs masked, or None if no CDRs found
        """
        if cdr_names is None:
            # Extract all CDRs and mask them
            cdrs = self.cdr_processor.extract_cdr_sequences(sequence)
            if not cdrs:
                logger.warning("No CDRs found in sequence")
                return None
            
            cdr_names = list(cdrs.keys())
        
        masked_sequence = sequence
        for cdr_name in cdr_names:
            masked_sequence = self.cdr_processor.mask_cdr_in_sequence(masked_sequence, cdr_name)
            if masked_sequence is None:
                return None
                
        return masked_sequence
    
    def predict_masked_cdrs(self, masked_sequence: str, cdr_names: List[str] = None, num_predictions: int = 5) -> Dict[str, List[str]]:
        """
        Predict masked CDR sequences in an antibody sequence.
        
        Args:
            masked_sequence: Antibody sequence with masked CDRs
            cdr_names: List of CDR names that were masked
            num_predictions: Number of predictions to generate per CDR
            
        Returns:
            Dictionary mapping CDR names to lists of predicted sequences
        """
        if "<mask>" not in masked_sequence:
            logger.warning("No mask tokens found in sequence")
            return {}
        
        if cdr_names is None:
            # Try to infer which CDRs are masked
            cdr_names = []
            cdrs = self.cdr_processor.extract_cdr_sequences(masked_sequence.replace("<mask>", "X" * 10))
            for cdr in cdrs:
                if "X" in cdrs[cdr]:
                    cdr_names.append(cdr)
        
        if not cdr_names:
            logger.warning("No CDR names provided or detected")
            return {}
        
        # Replace mask token with ESM-2 mask token
        masked_sequence = masked_sequence.replace("<mask>", "<mask>")
        
        # Tokenize sequence
        inputs = self.tokenizer(masked_sequence, return_tensors="pt")
        inputs = {k: v.to(self.device) for k, v in inputs.items()}
        
        # Get mask token positions
        mask_positions = (inputs["input_ids"] == self.tokenizer.mask_token_id).nonzero(as_tuple=True)[1]
        
        if len(mask_positions) == 0:
            logger.warning("No mask tokens found after tokenization")
            return {}
        
        # Forward pass through the model
        with torch.no_grad():
            outputs = self.model(**inputs)
            logits = outputs.logits
        
        # For each masked position, get top-k predictions
        all_predictions = {}
        
        # Group predictions by CDR
        for i, cdr_name in enumerate(cdr_names):
            if i >= len(mask_positions):
                logger.warning(f"Not enough mask positions for CDR {cdr_name}")
                continue
                
            position = mask_positions[i]
            token_logits = logits[0, position, :]
            
            # Get top-k predictions
            topk_probs, topk_indices = torch.topk(F.softmax(token_logits, dim=0), num_predictions)
            
            # Convert token IDs to amino acids
            predictions = []
            for j in range(num_predictions):
                token_id = topk_indices[j].item()
                amino_acid = self.tokenizer.decode([token_id]).strip()
                predictions.append(amino_acid)
                
            all_predictions[cdr_name] = predictions
        
        return all_predictions
    
    def generate_antibody(self, template_sequence: str, cdrs_to_modify: List[str], num_variants: int = 5) -> List[str]:
        """
        Generate antibody variants by modifying specific CDRs.
        
        Args:
            template_sequence: Template antibody sequence
            cdrs_to_modify: List of CDR names to modify
            num_variants: Number of variants to generate
            
        Returns:
            List of generated antibody sequences
        """
        if not cdrs_to_modify:
            logger.warning("No CDRs specified for modification")
            return []
        
        # Extract CDRs from template sequence
        cdrs = self.cdr_processor.extract_cdr_sequences(template_sequence)
        if not any(cdr in cdrs for cdr in cdrs_to_modify):
            logger.warning("None of the specified CDRs found in template sequence")
            return []
        
        variants = []
        attempts = 0
        max_attempts = self.config.antibody_generation["max_generation_attempts"]
        
        while len(variants) < num_variants and attempts < max_attempts:
            attempts += 1
            
            # Create a copy of the template sequence
            variant_sequence = template_sequence
            
            # For each CDR to modify, mask it and predict new sequences
            for cdr_name in cdrs_to_modify:
                if cdr_name not in cdrs:
                    continue
                    
                # Mask the CDR
                masked_sequence = self.cdr_processor.mask_cdr_in_sequence(variant_sequence, cdr_name)
                if masked_sequence is None:
                    continue
                    
                # Predict new CDR sequence
                predictions = self.predict_masked_cdrs(masked_sequence, [cdr_name], num_predictions=3)
                if not predictions or cdr_name not in predictions:
                    continue
                    
                # Select a prediction (based on diversity settings)
                if attempts % 3 == 0:  # Every third attempt, use more randomness
                    prediction_idx = random.choice(range(min(3, len(predictions[cdr_name]))))
                else:
                    prediction_idx = 0  # Use the top prediction
                    
                new_cdr = predictions[cdr_name][prediction_idx]
                
                # Replace the CDR with the new sequence
                variant_sequence = self.cdr_processor.replace_cdr_in_sequence(variant_sequence, cdr_name, new_cdr)
                if variant_sequence is None:
                    break
            
            # If we successfully created a variant, add it to the list
            if variant_sequence and variant_sequence != template_sequence:
                # Check if this variant is unique
                if variant_sequence not in variants:
                    variants.append(variant_sequence)
        
        logger.info(f"Generated {len(variants)} antibody variants after {attempts} attempts")
        return variants
    
    def predict_binding_score(self, toxin_sequence: str, antibody_sequence: str) -> float:
        """
        Predict binding score between toxin and antibody sequences.
        
        Args:
            toxin_sequence: Toxin amino acid sequence
            antibody_sequence: Antibody amino acid sequence
            
        Returns:
            Predicted binding score (0-1)
        """
        if not hasattr(self, "query_proj") or not hasattr(self, "binding_layers"):
            logger.warning("Binding prediction head not initialized")
            # Fallback to a simpler method
            return self._naive_binding_prediction(toxin_sequence, antibody_sequence)
        
        # Encode sequences
        toxin_embedding = self.encode_sequence(toxin_sequence)
        antibody_embedding = self.encode_sequence(antibody_sequence)
        
        # Cross-attention mechanism
        toxin_query = self.query_proj(toxin_embedding)
        antibody_key = self.key_proj(antibody_embedding)
        antibody_value = self.value_proj(antibody_embedding)
        
        # Compute attention scores
        attention_scores = torch.matmul(toxin_query, antibody_key.transpose(-1, -2)) / np.sqrt(self.config.model["hidden_dim"])
        attention_weights = F.softmax(attention_scores, dim=-1)
        
        # Compute context vector
        context_vector = torch.matmul(attention_weights, antibody_value)
        
        # Concatenate toxin query and context vector
        combined = torch.cat([toxin_query, context_vector], dim=1)
        
        # Predict binding score
        binding_score = self.binding_layers(combined)
        
        return binding_score.item()
    
    def _naive_binding_prediction(self, toxin_sequence: str, antibody_sequence: str) -> float:
        """
        Simple binding prediction based on embeddings similarity.
        
        Args:
            toxin_sequence: Toxin amino acid sequence
            antibody_sequence: Antibody amino acid sequence
            
        Returns:
            Binding score (0-1)
        """
        # Encode sequences
        toxin_embedding = self.encode_sequence(toxin_sequence)
        antibody_embedding = self.encode_sequence(antibody_sequence)
        
        # Extract CDRs
        cdrs = self.cdr_processor.extract_cdr_sequences(antibody_sequence)
        
        # If we found CDRs, use weighted similarity with higher weight on CDR regions
        if cdrs:
            # Create a new antibody sequence with only CDRs
            cdr_only_sequence = ''.join(cdrs.values())
            if len(cdr_only_sequence) > 0:
                cdr_embedding = self.encode_sequence(cdr_only_sequence)
                
                # Compute weighted similarity (70% CDR, 30% full antibody)
                cdr_sim = F.cosine_similarity(toxin_embedding, cdr_embedding)
                full_sim = F.cosine_similarity(toxin_embedding, antibody_embedding)
                similarity = 0.7 * cdr_sim + 0.3 * full_sim
                
                # Scale to 0-1 range
                binding_score = (similarity + 1) / 2
                return binding_score.item()
        
        # If no CDRs found, use cosine similarity between full sequences
        similarity = F.cosine_similarity(toxin_embedding, antibody_embedding)
        
        # Scale to 0-1 range
        binding_score = (similarity + 1) / 2
        return binding_score.item()
    
    def compute_sequence_score(self, sequence: str) -> float:
        """
        Compute a quality/naturalness score for an antibody sequence.
        
        Args:
            sequence: Antibody sequence
            
        Returns:
            Sequence score (0-1)
        """
        # Clean sequence
        sequence = re.sub(r"[^A-Za-z]", "", sequence).upper()
        
        # Tokenize sequence
        inputs = self.tokenizer(sequence, return_tensors="pt")
        inputs = {k: v.to(self.device) for k, v in inputs.items()}
        
        # Calculate perplexity (measure of how "natural" the sequence is)
        with torch.no_grad():
            outputs = self.model(**inputs)
            logits = outputs.logits
            
            # Calculate perplexity (excluding special tokens)
            shift_logits = logits[0, :-1, :].contiguous()
            shift_labels = inputs["input_ids"][0, 1:].contiguous()
            
            loss_fn = torch.nn.CrossEntropyLoss(reduction='mean')
            loss = loss_fn(shift_logits, shift_labels)
            perplexity = torch.exp(loss).item()
            
            # Convert perplexity to a score between 0 and 1 (lower perplexity = higher score)
            score = 1.0 / (1.0 + perplexity / 10.0)
            
            return score
    
    def save_checkpoint(self, save_path: str) -> None:
        """
        Save model checkpoint.
        
        Args:
            save_path: Path to save the checkpoint
        """
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        
        # Save binding prediction head if available
        binding_head_state = None
        if hasattr(self, "query_proj") and hasattr(self, "binding_layers"):
            binding_head_state = {
                "query_proj": self.query_proj.state_dict(),
                "key_proj": self.key_proj.state_dict(),
                "value_proj": self.value_proj.state_dict(),
                "binding_layers": self.binding_layers.state_dict()
            }
        
        # Save configuration and model states
        checkpoint = {
            "config": self.config.to_dict(),
            "model_state": self.model.state_dict(),
            "base_model_state": self.base_model.state_dict(),
            "binding_head_state": binding_head_state
        }
        
        torch.save(checkpoint, save_path)
        logger.info(f"Model checkpoint saved to {save_path}")
    
    def load_checkpoint(self, checkpoint_path: str) -> None:
        """
        Load model from checkpoint.
        
        Args:
            checkpoint_path: Path to the checkpoint
        """
        if not os.path.exists(checkpoint_path):
            logger.error(f"Checkpoint file not found: {checkpoint_path}")
            raise FileNotFoundError(f"Checkpoint file not found: {checkpoint_path}")
            
        try:
            checkpoint = torch.load(checkpoint_path, map_location=self.device)
            
            # Load configuration
            if "config" in checkpoint:
                self.config = ESM2Config.from_dict(checkpoint["config"])
            
            # Load model states
            if "model_state" in checkpoint:
                self.model.load_state_dict(checkpoint["model_state"])
            
            if "base_model_state" in checkpoint:
                self.base_model.load_state_dict(checkpoint["base_model_state"])
            
            # Load binding prediction head if available
            if "binding_head_state" in checkpoint and checkpoint["binding_head_state"] is not None:
                if not hasattr(self, "query_proj") or not hasattr(self, "binding_layers"):
                    self._init_binding_prediction_head()
                
                self.query_proj.load_state_dict(checkpoint["binding_head_state"]["query_proj"])
                self.key_proj.load_state_dict(checkpoint["binding_head_state"]["key_proj"])
                self.value_proj.load_state_dict(checkpoint["binding_head_state"]["value_proj"])
                self.binding_layers.load_state_dict(checkpoint["binding_head_state"]["binding_layers"])
            
            logger.info(f"Model loaded successfully from {checkpoint_path}")
        except Exception as e:
            logger.error(f"Error loading checkpoint: {e}")
            raise RuntimeError(f"Failed to load checkpoint: {e}")