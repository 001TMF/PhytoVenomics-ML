#!/usr/bin/env python
# utils/model_utils.py

import logging
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Dict, List, Any, Optional, Union, Tuple
import os
from pathlib import Path
import pickle
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
import joblib

logger = logging.getLogger("Phytovenomics.ModelUtils")

class ProteinLanguageModel(nn.Module):
    """Protein language model for encoding and generating protein sequences"""
    
    def __init__(self, embed_dim=512, num_heads=8, num_layers=6, dropout=0.1):
        super(ProteinLanguageModel, self).__init__()
        self.embed_dim = embed_dim
        self.vocab_size = 21  # 20 amino acids + padding
        
        # Embedding layer
        self.embedding = nn.Embedding(self.vocab_size, embed_dim)
        self.pos_encoding = nn.Parameter(torch.zeros(1, 1024, embed_dim))
        
        # Transformer layers
        encoder_layers = nn.TransformerEncoderLayer(
            d_model=embed_dim, nhead=num_heads, dropout=dropout, batch_first=True
        )
        self.transformer = nn.TransformerEncoder(encoder_layers, num_layers=num_layers)
        
        # Output projection
        self.fc_out = nn.Linear(embed_dim, self.vocab_size)
        
        logger.info(f"Initialized ProteinLanguageModel with {embed_dim} dimensions")
    
    def forward(self, x, mask=None):
        """Forward pass through the model"""
        seq_len = x.size(1)
        x = self.embedding(x) + self.pos_encoding[:, :seq_len, :]
        
        if mask is not None:
            x = self.transformer(x, src_key_padding_mask=mask)
        else:
            x = self.transformer(x)
            
        x = self.fc_out(x)
        return x
    
    def encode_sequence(self, seq_tensor):
        """Encode a protein sequence to embedding space"""
        with torch.no_grad():
            seq_embed = self.embedding(seq_tensor)
            seq_len = seq_tensor.size(1)
            seq_embed = seq_embed + self.pos_encoding[:, :seq_len, :]
            encoded = self.transformer(seq_embed)
            # Return mean embedding across sequence length
            return encoded.mean(dim=1)


class AntibodyBindingPredictor(nn.Module):
    """Neural network for predicting antibody-toxin binding affinity"""
    
    def __init__(self, embed_dim=512, hidden_dims=[512, 256, 128]):
        super(AntibodyBindingPredictor, self).__init__()
        
        self.embed_dim = embed_dim
        layers = []
        
        input_dim = embed_dim * 2  # Concatenated antibody and toxin embeddings
        
        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(input_dim, hidden_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(0.2))
            input_dim = hidden_dim
            
        layers.append(nn.Linear(hidden_dims[-1], 1))  # Output binding energy
        
        self.network = nn.Sequential(*layers)
        logger.info(f"Initialized AntibodyBindingPredictor with {len(hidden_dims)} hidden layers")
    
    def forward(self, antibody_embed, toxin_embed):
        """Forward pass to predict binding affinity
        
        Args:
            antibody_embed: Tensor of antibody embeddings [batch_size, embed_dim]
            toxin_embed: Tensor of toxin embeddings [batch_size, embed_dim]
            
        Returns:
            Predicted binding energy (lower = stronger binding)
        """
        combined = torch.cat([antibody_embed, toxin_embed], dim=1)
        binding_energy = self.network(combined)
        return binding_energy


class EpitopePredictor(nn.Module):
    """Model for predicting epitopes on toxin sequences"""
    
    def __init__(self, embed_dim=512, hidden_dim=256, num_layers=2):
        super(EpitopePredictor, self).__init__()
        
        self.embed_dim = embed_dim
        
        # Bidirectional LSTM for sequence context
        self.lstm = nn.LSTM(
            input_size=embed_dim,
            hidden_size=hidden_dim,
            num_layers=num_layers,
            batch_first=True,
            bidirectional=True
        )
        
        # Output layer
        self.classifier = nn.Linear(hidden_dim * 2, 1)  # *2 for bidirectional
        self.sigmoid = nn.Sigmoid()
        
        logger.info(f"Initialized EpitopePredictor with LSTM hidden dim {hidden_dim}")
    
    def forward(self, x):
        """Forward pass to predict epitope probability for each position"""
        lstm_out, _ = self.lstm(x)
        logits = self.classifier(lstm_out)
        probs = self.sigmoid(logits)
        return probs.squeeze(-1)


class ModelUtils:
    """Utility functions for ML models in the Phytovenomics platform."""
    
    @staticmethod
    def check_cuda_availability() -> bool:
        """Check if CUDA is available for GPU acceleration.
        
        Returns:
            Boolean indicating CUDA availability
        """
        cuda_available = torch.cuda.is_available()
        cuda_device_count = torch.cuda.device_count() if cuda_available else 0
        
        logger.info(f"CUDA available: {cuda_available}, Device count: {cuda_device_count}")
        
        if cuda_available:
            for i in range(cuda_device_count):
                device_name = torch.cuda.get_device_name(i)
                logger.info(f"  Device {i}: {device_name}")
                
        return cuda_available
    
    @staticmethod
    def set_random_seed(seed: int = 42):
        """Set random seed for reproducibility.
        
        Args:
            seed: Random seed value
        """
        logger.info(f"Setting random seed to {seed}")
        np.random.seed(seed)
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    
    @staticmethod
    def encode_protein_sequence(sequence: str, max_length: int = 512) -> np.ndarray:
        """Encode a protein sequence using one-hot encoding.
        
        Args:
            sequence: Protein sequence as a string
            max_length: Maximum sequence length
            
        Returns:
            One-hot encoded sequence as numpy array
        """
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}
        
        # Truncate or pad sequence
        if len(sequence) > max_length:
            sequence = sequence[:max_length]
        else:
            sequence = sequence + "_" * (max_length - len(sequence))
        
        # Create one-hot encoding
        encoding = np.zeros((max_length, len(amino_acids)))
        for i, aa in enumerate(sequence):
            if aa in aa_to_idx:  # Skip padding character
                encoding[i, aa_to_idx[aa]] = 1.0
        
        return encoding
    
    @staticmethod
    def encode_amino_acids(sequence: str, aa_to_idx: Dict[str, int]) -> torch.Tensor:
        """Convert an amino acid sequence to indices for model input.
        
        Args:
            sequence: Amino acid sequence
            aa_to_idx: Mapping from amino acids to indices
            
        Returns:
            Tensor of amino acid indices
        """
        indices = [aa_to_idx.get(aa.upper(), 0) for aa in sequence]
        return torch.tensor(indices, dtype=torch.long)
    
    @staticmethod
    def sliding_window(sequence: str, window_size: int = 21) -> List[str]:
        """Generate sliding windows from a sequence.
        
        Args:
            sequence: Amino acid sequence
            window_size: Size of the window
            
        Returns:
            List of window sequences
        """
        windows = []
        half_size = window_size // 2
        
        # Pad sequence
        padded_sequence = "X" * half_size + sequence + "X" * half_size
        
        # Generate windows
        for i in range(len(sequence)):
            window = padded_sequence[i:i+window_size]
            windows.append(window)
            
        return windows
    
    @staticmethod
    def calculate_sequence_similarity(seq1: str, seq2: str) -> float:
        """Calculate sequence similarity between two protein sequences.
        
        Args:
            seq1: First protein sequence
            seq2: Second protein sequence
            
        Returns:
            Similarity score (0-1)
        """
        # Ensure equal length for comparison
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]
        
        # Count matching positions
        matches = sum(a == b for a, b in zip(seq1, seq2))
        
        # Calculate similarity
        return matches / min_len if min_len > 0 else 0
    
    @staticmethod
    def batch_sequences(sequences: List[str], batch_size: int = 32) -> List[List[str]]:
        """Batch sequences for more efficient processing.
        
        Args:
            sequences: List of protein sequences
            batch_size: Size of each batch
            
        Returns:
            List of sequence batches
        """
        return [sequences[i:i + batch_size] for i in range(0, len(sequences), batch_size)]
    
    @staticmethod
    def load_model_weights(model, weights_path: str) -> bool:
        """Load pre-trained model weights.
        
        Args:
            model: Model object
            weights_path: Path to weights file
            
        Returns:
            Boolean indicating success
        """
        try:
            if not os.path.exists(weights_path):
                logger.error(f"Weights file not found: {weights_path}")
                return False
            
            if torch.cuda.is_available():
                state_dict = torch.load(weights_path)
            else:
                state_dict = torch.load(weights_path, map_location=torch.device('cpu'))
                
            model.load_state_dict(state_dict)
            logger.info(f"Model weights loaded from {weights_path}")
            return True
        except Exception as e:
            logger.error(f"Error loading model weights: {e}")
            return False
    
    @staticmethod
    def save_model_weights(model, save_path: str):
        """Save model weights to file.
        
        Args:
            model: Model object
            save_path: Path to save weights file
        """
        try:
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            torch.save(model.state_dict(), save_path)
            logger.info(f"Model weights saved to {save_path}")
        except Exception as e:
            logger.error(f"Error saving model weights: {e}")
            
    @staticmethod
    def get_device() -> torch.device:
        """Get device to use for model operations.
        
        Returns:
            Torch device object
        """
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    @staticmethod
    def calculate_model_memory_usage(model, input_size: Tuple) -> float:
        """Estimate memory usage of a model in GB.
        
        Args:
            model: PyTorch model
            input_size: Size of input tensor
            
        Returns:
            Estimated memory usage in GB
        """
        try:
            # Calculate parameter memory
            param_size = sum(p.numel() for p in model.parameters()) * 4  # bytes for float32
            
            # Calculate input and output memory
            input_size_bytes = np.prod(input_size) * 4  # bytes for float32
            
            # Estimate total memory usage
            total_gb = (param_size + input_size_bytes * 3) / (1024**3)  # Convert to GB
            
            logger.info(f"Estimated model memory usage: {total_gb:.2f} GB")
            return total_gb
        except Exception as e:
            logger.error(f"Error calculating model memory usage: {e}")
            return 0.0
    
    @staticmethod
    def train_binding_predictor(model, train_loader, val_loader, optimizer, num_epochs=10, device='cpu'):
        """Train the binding affinity prediction model.
        
        Args:
            model: AntibodyBindingPredictor model
            train_loader: DataLoader for training data
            val_loader: DataLoader for validation data
            optimizer: Optimizer instance
            num_epochs: Number of training epochs
            device: Device to train on ('cpu' or 'cuda')
            
        Returns:
            Dictionary with training metrics
        """
        logger.info(f"Training binding predictor model for {num_epochs} epochs on {device}")
        model.to(device)
        best_val_loss = float('inf')
        metrics = {'train_loss': [], 'val_loss': []}
        
        for epoch in range(num_epochs):
            # Training phase
            model.train()
            train_loss = 0.0
            
            for ab_embeds, toxin_embeds, affinities in train_loader:
                ab_embeds = ab_embeds.to(device)
                toxin_embeds = toxin_embeds.to(device)
                affinities = affinities.to(device)
                
                optimizer.zero_grad()
                predictions = model(ab_embeds, toxin_embeds)
                loss = F.mse_loss(predictions, affinities)
                loss.backward()
                optimizer.step()
                
                train_loss += loss.item() * ab_embeds.size(0)
            
            train_loss /= len(train_loader.dataset)
            metrics['train_loss'].append(train_loss)
            
            # Validation phase
            model.eval()
            val_loss = 0.0
            
            with torch.no_grad():
                for ab_embeds, toxin_embeds, affinities in val_loader:
                    ab_embeds = ab_embeds.to(device)
                    toxin_embeds = toxin_embeds.to(device)
                    affinities = affinities.to(device)
                    
                    predictions = model(ab_embeds, toxin_embeds)
                    loss = F.mse_loss(predictions, affinities)
                    
                    val_loss += loss.item() * ab_embeds.size(0)
                
                val_loss /= len(val_loader.dataset)
                metrics['val_loss'].append(val_loss)
                
                logger.info(f"Epoch {epoch+1}/{num_epochs} - Train loss: {train_loss:.4f}, Val loss: {val_loss:.4f}")
                
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    # Save best model checkpoint
                    torch.save(model.state_dict(), 'models/best_binding_predictor.pt')
                    logger.info(f"Saved best model with val_loss: {val_loss:.4f}")
        
        return metrics
    
    @staticmethod
    def train_epitope_predictor(model, train_loader, val_loader, optimizer, num_epochs=10, device='cpu'):
        """Train the epitope prediction model.
        
        Args:
            model: EpitopePredictor model
            train_loader: DataLoader for training data
            val_loader: DataLoader for validation data
            optimizer: Optimizer instance
            num_epochs: Number of training epochs
            device: Device to train on ('cpu' or 'cuda')
            
        Returns:
            Dictionary with training metrics
        """
        logger.info(f"Training epitope predictor model for {num_epochs} epochs on {device}")
        model.to(device)
        best_val_loss = float('inf')
        metrics = {'train_loss': [], 'val_loss': []}
        
        for epoch in range(num_epochs):
            # Training phase
            model.train()
            train_loss = 0.0
            
            for toxin_embeds, epitope_labels in train_loader:
                toxin_embeds = toxin_embeds.to(device)
                epitope_labels = epitope_labels.to(device)
                
                optimizer.zero_grad()
                predictions = model(toxin_embeds)
                loss = F.binary_cross_entropy(predictions, epitope_labels)
                loss.backward()
                optimizer.step()
                
                train_loss += loss.item() * toxin_embeds.size(0)
            
            train_loss /= len(train_loader.dataset)
            metrics['train_loss'].append(train_loss)
            
            # Validation phase
            model.eval()
            val_loss = 0.0
            
            with torch.no_grad():
                for toxin_embeds, epitope_labels in val_loader:
                    toxin_embeds = toxin_embeds.to(device)
                    epitope_labels = epitope_labels.to(device)
                    
                    predictions = model(toxin_embeds)
                    loss = F.binary_cross_entropy(predictions, epitope_labels)
                    
                    val_loss += loss.item() * toxin_embeds.size(0)
                
                val_loss /= len(val_loader.dataset)
                metrics['val_loss'].append(val_loss)
                
                logger.info(f"Epoch {epoch+1}/{num_epochs} - Train loss: {train_loss:.4f}, Val loss: {val_loss:.4f}")
                
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    # Save best model checkpoint
                    torch.save(model.state_dict(), 'models/best_epitope_predictor.pt')
                    logger.info(f"Saved best model with val_loss: {val_loss:.4f}")
        
        return metrics
    
    @staticmethod
    def train_cdr_generator(model, train_loader, val_loader, optimizer, num_epochs=10, device='cpu'):
        """Train the CDR generator model.
        
        Args:
            model: ProteinLanguageModel model
            train_loader: DataLoader for training data
            val_loader: DataLoader for validation data
            optimizer: Optimizer instance
            num_epochs: Number of training epochs
            device: Device to train on ('cpu' or 'cuda')
            
        Returns:
            Dictionary with training metrics
        """
        logger.info(f"Training CDR generator model for {num_epochs} epochs on {device}")
        model.to(device)
        best_val_loss = float('inf')
        metrics = {'train_loss': [], 'val_loss': []}
        
        for epoch in range(num_epochs):
            # Training phase
            model.train()
            train_loss = 0.0
            
            for src_seq, tgt_seq, mask in train_loader:
                src_seq = src_seq.to(device)
                tgt_seq = tgt_seq.to(device)
                mask = mask.to(device) if mask is not None else None
                
                optimizer.zero_grad()
                logits = model(src_seq, mask)
                
                # Reshape for cross entropy
                batch_size, seq_len, vocab_size = logits.size()
                loss = F.cross_entropy(
                    logits.view(batch_size * seq_len, vocab_size),
                    tgt_seq.view(-1),
                    ignore_index=0  # Ignore padding
                )
                
                loss.backward()
                optimizer.step()
                
                train_loss += loss.item() * batch_size
            
            train_loss /= len(train_loader.dataset)
            metrics['train_loss'].append(train_loss)
            
            # Validation phase
            model.eval()
            val_loss = 0.0
            
            with torch.no_grad():
                for src_seq, tgt_seq, mask in val_loader:
                    src_seq = src_seq.to(device)
                    tgt_seq = tgt_seq.to(device)
                    mask = mask.to(device) if mask is not None else None
                    
                    logits = model(src_seq, mask)
                    
                    # Reshape for cross entropy
                    batch_size, seq_len, vocab_size = logits.size()
                    loss = F.cross_entropy(
                        logits.view(batch_size * seq_len, vocab_size),
                        tgt_seq.view(-1),
                        ignore_index=0  # Ignore padding
                    )
                    
                    val_loss += loss.item() * batch_size
                
                val_loss /= len(val_loader.dataset)
                metrics['val_loss'].append(val_loss)
                
                logger.info(f"Epoch {epoch+1}/{num_epochs} - Train loss: {train_loss:.4f}, Val loss: {val_loss:.4f}")
                
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    # Save best model checkpoint
                    torch.save(model.state_dict(), 'models/best_cdr_generator.pt')
                    logger.info(f"Saved best model with val_loss: {val_loss:.4f}")
        
        return metrics
        
    @staticmethod
    def train_random_forest_binding_model(X_train, y_train, X_val=None, y_val=None, n_estimators=100):
        """Train a Random Forest model for binding affinity prediction.
        
        Args:
            X_train: Training features
            y_train: Training labels (binding affinities)
            X_val: Validation features (optional)
            y_val: Validation labels (optional)
            n_estimators: Number of trees in the forest
            
        Returns:
            Trained model and scaler
        """
        logger.info(f"Training Random Forest binding model with {n_estimators} trees")
        
        # Scale features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        
        # Train model
        model = RandomForestRegressor(
            n_estimators=n_estimators, 
            max_depth=10,
            min_samples_split=5,
            random_state=42,
            n_jobs=-1
        )
        
        model.fit(X_train_scaled, y_train)
        
        # Evaluate if validation data provided
        if X_val is not None and y_val is not None:
            X_val_scaled = scaler.transform(X_val)
            y_pred = model.predict(X_val_scaled)
            
            mse = mean_squared_error(y_val, y_pred)
            r2 = r2_score(y_val, y_pred)
            
            logger.info(f"Validation MSE: {mse:.4f}, R² Score: {r2:.4f}")
        
        # Save model and scaler
        os.makedirs('models', exist_ok=True)
        joblib.dump(model, 'models/rf_binding_model.joblib')
        joblib.dump(scaler, 'models/rf_binding_scaler.joblib')
        
        return model, scaler
    
    @staticmethod
    def create_epitope_prediction_datasets(toxins_df, epitopes_df, val_split=0.2, negative_ratio=1.0):
        """Create datasets for epitope prediction model training.
        
        Args:
            toxins_df: DataFrame with toxin sequences
            epitopes_df: DataFrame with known epitopes
            val_split: Fraction of data to use for validation
            negative_ratio: Ratio of negative to positive examples
            
        Returns:
            Dictionary with training and validation datasets
        """
        logger.info(f"Creating epitope prediction datasets from {len(toxins_df)} toxins and {len(epitopes_df)} epitopes")
        
        # Initialize lists for features and labels
        X_sequences = []
        y_labels = []
        
        # Process each toxin
        for _, toxin_row in toxins_df.iterrows():
            toxin_id = toxin_row['id']
            sequence = toxin_row['sequence']
            
            # Get epitopes for this toxin
            toxin_epitopes = epitopes_df[epitopes_df['toxin_id'] == toxin_id]
            
            # Create positive examples from known epitopes
            pos_positions = set()
            for _, epitope_row in toxin_epitopes.iterrows():
                start = int(epitope_row['start'])
                end = int(epitope_row['end'])
                # Mark all positions in this epitope
                for pos in range(start, end + 1):
                    if pos < len(sequence):
                        pos_positions.add(pos)
            
            # Create sliding window examples
            window_size = 21  # Odd number for central residue
            half_window = window_size // 2
            
            # Pad sequence for edge cases
            padded_seq = "X" * half_window + sequence + "X" * half_window
            
            # Generate windows for each position
            for i in range(len(sequence)):
                # Extract window centered at position i
                window = padded_seq[i:i + window_size]
                if len(window) == window_size:  # Ensure full window
                    X_sequences.append(window)
                    # Label as 1 if center position is in an epitope
                    y_labels.append(1.0 if i in pos_positions else 0.0)
        
        # Balance dataset if needed
        X_pos = [seq for seq, label in zip(X_sequences, y_labels) if label > 0.5]
        y_pos = [1.0] * len(X_pos)
        X_neg = [seq for seq, label in zip(X_sequences, y_labels) if label < 0.5]
        y_neg = [0.0] * len(X_neg)
        
        # Subsample negative examples if needed
        if negative_ratio < 1.0 and len(X_neg) > len(X_pos) * negative_ratio:
            neg_sample_size = int(len(X_pos) * negative_ratio)
            indices = np.random.choice(len(X_neg), neg_sample_size, replace=False)
            X_neg = [X_neg[i] for i in indices]
            y_neg = [0.0] * len(X_neg)
        
        # Combine positive and negative examples
        X_combined = X_pos + X_neg
        y_combined = y_pos + y_neg
        
        # Shuffle the data
        indices = np.arange(len(X_combined))
        np.random.shuffle(indices)
        X_shuffled = [X_combined[i] for i in indices]
        y_shuffled = [y_combined[i] for i in indices]
        
        # Split into training and validation
        split_idx = int(len(X_shuffled) * (1 - val_split))
        X_train, X_val = X_shuffled[:split_idx], X_shuffled[split_idx:]
        y_train, y_val = y_shuffled[:split_idx], y_shuffled[split_idx:]
        
        logger.info(f"Created dataset with {len(X_train)} training samples and {len(X_val)} validation samples")
        logger.info(f"Positive examples: {sum(y_train) + sum(y_val)}, Negative examples: {len(y_train) + len(y_val) - sum(y_train) - sum(y_val)}")
        
        return {
            'X_train': X_train, 
            'y_train': y_train, 
            'X_val': X_val, 
            'y_val': y_val
        }