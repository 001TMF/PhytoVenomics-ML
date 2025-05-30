#!/usr/bin/env python
# utils/training_utils.py

import os
import random
import logging
import numpy as np
import torch
from typing import Union, Optional

logger = logging.getLogger("Phytovenomics.TrainingUtils")

def set_seed(seed: int = 42) -> None:
    """
    Set random seed for reproducibility across all libraries.
    
    Args:
        seed: Random seed value
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    # When running on the CuDNN backend, two further options must be set
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    # Set a fixed value for the hash seed
    os.environ["PYTHONHASHSEED"] = str(seed)
    
    logger.info(f"Random seed set to {seed} for reproducibility")

def setup_device(device_str: str = "auto") -> torch.device:
    """
    Set up device for computation (CPU or GPU).
    
    Args:
        device_str: Device string ("auto", "cpu", "cuda", "cuda:0", etc.)
        
    Returns:
        PyTorch device
    """
    if device_str == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(device_str)
    
    if device.type == "cuda":
        logger.info(f"Using GPU: {torch.cuda.get_device_name(device.index or 0)}")
        logger.info(f"Available memory: {torch.cuda.get_device_properties(device.index or 0).total_memory / 1e9:.2f} GB")
    else:
        logger.info("Using CPU")
        
    return device

def manage_checkpoints(
    model: torch.nn.Module, 
    save_dir: str, 
    epoch: int, 
    loss: float, 
    optimizer: Optional[torch.optim.Optimizer] = None,
    scheduler: Optional[torch.optim.lr_scheduler._LRScheduler] = None,
    save_optimizer: bool = True,
    keep_latest_k: int = 5,
    best_only: bool = False,
    is_best: bool = False,
    model_name: str = "esm2_cdr"
) -> str:
    """
    Save model checkpoint and manage checkpoint directory.
    
    Args:
        model: PyTorch model to save
        save_dir: Directory to save checkpoints
        epoch: Current epoch number
        loss: Validation loss
        optimizer: Optimizer (optional)
        scheduler: Learning rate scheduler (optional)
        save_optimizer: Whether to save optimizer state
        keep_latest_k: Number of latest checkpoints to keep
        best_only: Whether to save only the best checkpoint
        is_best: Whether current model is the best so far
        model_name: Name prefix for checkpoint files
        
    Returns:
        Path to the saved checkpoint
    """
    os.makedirs(save_dir, exist_ok=True)
    
    # Prepare checkpoint
    checkpoint = {
        "epoch": epoch,
        "loss": loss,
        "model_state_dict": model.state_dict()
    }
    
    if save_optimizer and optimizer is not None:
        checkpoint["optimizer_state_dict"] = optimizer.state_dict()
        
    if scheduler is not None:
        checkpoint["scheduler_state_dict"] = scheduler.state_dict()
    
    # Save checkpoint
    if best_only:
        # Save only if it's the best model
        if is_best:
            checkpoint_path = os.path.join(save_dir, f"{model_name}_best.pt")
            torch.save(checkpoint, checkpoint_path)
            logger.info(f"Saved best model checkpoint to {checkpoint_path}")
    else:
        # Regular checkpoint saving
        checkpoint_path = os.path.join(save_dir, f"{model_name}_epoch_{epoch}.pt")
        torch.save(checkpoint, checkpoint_path)
        logger.info(f"Saved checkpoint at epoch {epoch} to {checkpoint_path}")
        
        # Save best model separately if specified
        if is_best:
            best_path = os.path.join(save_dir, f"{model_name}_best.pt")
            torch.save(checkpoint, best_path)
            logger.info(f"Saved best model checkpoint to {best_path}")
            
        # Remove old checkpoints if needed
        if keep_latest_k > 0:
            checkpoints = [f for f in os.listdir(save_dir) if f.startswith(model_name) and f.endswith(".pt") and "best" not in f]
            checkpoints.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]) if x.split('_')[-1].split('.')[0].isdigit() else 0)
            
            for old_ckpt in checkpoints[:-keep_latest_k]:
                old_path = os.path.join(save_dir, old_ckpt)
                os.remove(old_path)
                logger.info(f"Removed old checkpoint: {old_path}")
    
    return checkpoint_path

def load_checkpoint(
    model: torch.nn.Module, 
    checkpoint_path: str, 
    optimizer: Optional[torch.optim.Optimizer] = None,
    scheduler: Optional[torch.optim.lr_scheduler._LRScheduler] = None,
    device: Optional[torch.device] = None
) -> dict:
    """
    Load model from checkpoint.
    
    Args:
        model: PyTorch model to load weights into
        checkpoint_path: Path to the checkpoint file
        optimizer: Optimizer to load state (optional)
        scheduler: Learning rate scheduler to load state (optional)
        device: Device to load the model to (optional)
        
    Returns:
        Dictionary with loaded checkpoint information
    """
    if not os.path.exists(checkpoint_path):
        logger.error(f"Checkpoint file not found: {checkpoint_path}")
        raise FileNotFoundError(f"Checkpoint file not found: {checkpoint_path}")
    
    # Load checkpoint with correct device mapping
    if device is None:
        device = next(model.parameters()).device
        
    checkpoint = torch.load(checkpoint_path, map_location=device)
    
    # Load model weights
    model.load_state_dict(checkpoint["model_state_dict"])
    
    # Load optimizer state if provided
    if optimizer is not None and "optimizer_state_dict" in checkpoint:
        optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
    
    # Load scheduler state if provided
    if scheduler is not None and "scheduler_state_dict" in checkpoint:
        scheduler.load_state_dict(checkpoint["scheduler_state_dict"])
    
    logger.info(f"Loaded checkpoint from {checkpoint_path} (epoch: {checkpoint.get('epoch', 'unknown')})")
    return checkpoint

def get_linear_schedule_with_warmup(
    optimizer: torch.optim.Optimizer, 
    num_warmup_steps: int, 
    num_training_steps: int
) -> torch.optim.lr_scheduler.LambdaLR:
    """
    Create a learning rate scheduler with warmup.
    
    Args:
        optimizer: Optimizer
        num_warmup_steps: Number of warmup steps
        num_training_steps: Total number of training steps
        
    Returns:
        Learning rate scheduler
    """
    def lr_lambda(current_step: int):
        if current_step < num_warmup_steps:
            # Linear warmup
            return float(current_step) / float(max(1, num_warmup_steps))
        # Linear decay after warmup
        return max(
            0.0, float(num_training_steps - current_step) / float(max(1, num_training_steps - num_warmup_steps))
        )
        
    return torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)

def mixed_precision_setup(use_fp16: bool = True) -> Optional[torch.cuda.amp.GradScaler]:
    """
    Set up mixed precision training.
    
    Args:
        use_fp16: Whether to use FP16 mixed precision
        
    Returns:
        Gradient scaler for mixed precision, or None if not using FP16
    """
    if use_fp16 and torch.cuda.is_available():
        scaler = torch.cuda.amp.GradScaler()
        logger.info("Mixed precision training enabled (FP16)")
        return scaler
    
    logger.info("Using full precision training (FP32)")
    return None

def get_model_size(model: torch.nn.Module) -> int:
    """
    Calculate model size in parameters.
    
    Args:
        model: PyTorch model
        
    Returns:
        Number of parameters in the model
    """
    return sum(p.numel() for p in model.parameters())

def log_gpu_memory_usage(device: torch.device) -> None:
    """
    Log GPU memory usage statistics.
    
    Args:
        device: PyTorch device to check
    """
    if device.type == "cuda":
        # Memory allocated in bytes
        allocated = torch.cuda.memory_allocated(device)
        # Maximum memory allocated in bytes
        max_allocated = torch.cuda.max_memory_allocated(device)
        # Memory reserved by the caching allocator in bytes
        reserved = torch.cuda.memory_reserved(device)
        
        logger.info(f"GPU memory usage - Allocated: {allocated/1e9:.2f}GB, "
                    f"Max allocated: {max_allocated/1e9:.2f}GB, "
                    f"Reserved: {reserved/1e9:.2f}GB")