"""
IgLM-based antibody sequence optimization

Based on: germinal-main/colabdesign/colabdesign/iglm/model.py
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Optional, List, Dict, Tuple

try:
    from iglm import IgLM
except ImportError:
    raise ImportError(
        "IgLM not installed. Install with: pip install iglm"
    )


class IgLMOptimizer(nn.Module):
    """
    IgLM-based antibody sequence optimizer.

    This class wraps the IgLM language model to evaluate and optimize
    antibody sequences for human-likeness and naturalness.

    Uses straight-through estimator for gradient computation.

    Based on: germinal-main/colabdesign/colabdesign/iglm/model.py
    """

    def __init__(
        self,
        model_name: str = "IgLM",
        chain_token: str = "[HEAVY]",
        species: str = "[HUMAN]",
        temperature: float = 1.0,
        cdr_positions: Optional[List[Tuple[int, int]]] = None,
        device: Optional[torch.device] = None,
        seed: int = 0,
    ):
        """
        Initialize IgLM optimizer.

        Args:
            model_name: IgLM model name ("IgLM" or "IgLM-S")
            chain_token: Chain token ("[HEAVY]" or "[LIGHT]")
            species: Species token ("[HUMAN]", "[MOUSE]", etc.)
            temperature: Softmax temperature for sequence sampling
            cdr_positions: List of (start, end) positions for CDR regions
            device: Torch device
            seed: Random seed
        """
        super().__init__()

        # Initialize IgLM
        self.iglm = IgLM(model_name=model_name)

        if device is not None:
            self.device = device
        else:
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        self.iglm.model.to(self.device)

        # Freeze IgLM parameters
        for param in self.iglm.model.parameters():
            param.requires_grad = False

        self.chain_token = chain_token
        self.species_token = species
        self.cdr_positions = cdr_positions
        self.temperature = temperature

        # Amino acid vocabulary
        self.amino_acids = [
            'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
        ]

        # Get amino acid token IDs
        aa_ids = []
        for aa in self.amino_acids:
            tid = self.iglm.tokenizer.convert_tokens_to_ids(aa)
            if tid == self.iglm.tokenizer.unk_token_id:
                raise ValueError(f"Unrecognized amino acid token: {aa}")
            aa_ids.append(tid)
        self.amino_acid_ids = torch.tensor(aa_ids, device=self.device)

        # Get special token IDs
        self.chain_id = self.iglm.tokenizer.convert_tokens_to_ids(chain_token)
        self.species_id = self.iglm.tokenizer.convert_tokens_to_ids(species)
        self.suffix_id = self.iglm.tokenizer.sep_token_id

        # Set seed
        if seed is not None:
            torch.manual_seed(seed)

    def forward(
        self,
        seq_logits: torch.Tensor,
        chain: str = "[HEAVY]"
    ) -> Tuple[torch.Tensor, float]:
        """
        Compute IgLM loss and log-likelihood for sequence.

        Args:
            seq_logits: Sequence logits (L, 20)
            chain: Chain type ("[HEAVY]" or "[LIGHT]")

        Returns:
            Tuple of (loss, log_likelihood)
        """
        # Apply softmax with temperature
        soft_probs = F.softmax(seq_logits / self.temperature, dim=-1)

        # Get hard indices
        hard_indices = soft_probs.argmax(dim=-1)
        hard_one_hot = F.one_hot(hard_indices, num_classes=soft_probs.size(-1)).float()

        # Straight-through estimator: forward uses hard, backward uses soft
        ste_probs = hard_one_hot + (soft_probs - soft_probs.detach())

        # Get amino acid embeddings
        embed_layer = self.iglm.model.get_input_embeddings()
        amino_embeds = embed_layer(self.amino_acid_ids)  # (20, embed_dim)

        # Compute variable region embeddings
        var_embeds = ste_probs @ amino_embeds  # (L, embed_dim)

        # Construct full sequence with prefix and suffix
        chain_id = self.iglm.tokenizer.convert_tokens_to_ids(chain)
        prefix_ids = torch.tensor([chain_id, self.species_id], device=self.device)
        prefix_embeds = embed_layer(prefix_ids)  # (2, embed_dim)

        suffix_ids = torch.tensor([self.suffix_id], device=self.device)
        suffix_embeds = embed_layer(suffix_ids)  # (1, embed_dim)

        # Concatenate: [prefix, variable, suffix]
        full_embeds = torch.cat([prefix_embeds, var_embeds, suffix_embeds], dim=0).unsqueeze(0)

        # Forward through IgLM
        outputs = self.iglm.model(inputs_embeds=full_embeds)
        logits = outputs.logits  # (1, total_length, vocab_size)

        # Construct target token IDs
        var_token_ids = self.amino_acid_ids[hard_indices]
        full_target_ids = torch.cat([prefix_ids, var_token_ids, suffix_ids], dim=0).unsqueeze(0)

        # Compute cross-entropy loss
        shift_logits = logits[:, :-1, :]
        shift_labels = full_target_ids[:, 1:]
        ce_loss = F.cross_entropy(
            shift_logits.reshape(-1, shift_logits.size(-1)),
            shift_labels.reshape(-1),
            reduction='none'
        )

        # Reshape and compute mean loss
        position_losses = ce_loss.reshape(shift_labels.shape)
        position_losses = position_losses[:, 1:-1]  # Exclude prefix/suffix

        total_loss = position_losses.mean()
        log_likelihood = -total_loss.item()

        return total_loss, log_likelihood

    def get_gradient(
        self,
        seq_logits: torch.Tensor,
        chain: str = "[HEAVY]"
    ) -> Tuple[torch.Tensor, float]:
        """
        Compute gradient of IgLM loss w.r.t. sequence logits.

        Args:
            seq_logits: Sequence logits (L, 20), requires_grad=True
            chain: Chain type

        Returns:
            Tuple of (gradient, log_likelihood)
        """
        ce_loss, ll = self.forward(seq_logits, chain=chain)

        # Compute gradient
        gradient = torch.autograd.grad(ce_loss, seq_logits)[0]

        return gradient.detach(), ll

    def score_sequence(
        self,
        sequence: str,
        chain: str = "[HEAVY]"
    ) -> Dict[str, float]:
        """
        Score a sequence using IgLM.

        Args:
            sequence: Amino acid sequence
            chain: Chain type

        Returns:
            Dictionary with log-likelihood and perplexity
        """
        # Convert sequence to logits (one-hot encoding)
        seq_indices = [self.amino_acids.index(aa) for aa in sequence]
        seq_logits = torch.zeros(len(sequence), 20, device=self.device)
        for i, idx in enumerate(seq_indices):
            seq_logits[i, idx] = 10.0  # High value for one-hot

        seq_logits.requires_grad = True

        # Compute log-likelihood
        with torch.no_grad():
            _, log_likelihood = self.forward(seq_logits, chain=chain)

        perplexity = np.exp(-log_likelihood / len(sequence))

        return {
            'log_likelihood': log_likelihood,
            'perplexity': perplexity,
            'avg_ll_per_residue': log_likelihood / len(sequence)
        }

    def score_cdr_regions(
        self,
        sequence: str,
        cdr_positions: List[Tuple[int, int]],
        chain: str = "[HEAVY]"
    ) -> Dict[str, Dict[str, float]]:
        """
        Score individual CDR regions.

        Args:
            sequence: Full antibody sequence
            cdr_positions: List of (start, end) positions for each CDR
            chain: Chain type

        Returns:
            Dictionary mapping CDR names to scores
        """
        cdr_names = ['H1', 'H2', 'H3'] if chain == "[HEAVY]" else ['L1', 'L2', 'L3']

        scores = {}
        for cdr_name, (start, end) in zip(cdr_names, cdr_positions):
            cdr_seq = sequence[start:end]
            cdr_score = self.score_sequence(cdr_seq, chain=chain)
            scores[cdr_name] = cdr_score

        return scores

    def suggest_mutations(
        self,
        sequence: str,
        num_suggestions: int = 5,
        chain: str = "[HEAVY]"
    ) -> List[Dict[str, any]]:
        """
        Suggest mutations to improve human-likeness.

        Args:
            sequence: Current sequence
            num_suggestions: Number of mutations to suggest
            chain: Chain type

        Returns:
            List of mutation suggestions with scores
        """
        # Convert to logits
        seq_indices = [self.amino_acids.index(aa) for aa in sequence]
        seq_logits = torch.zeros(len(sequence), 20, device=self.device, requires_grad=True)
        for i, idx in enumerate(seq_indices):
            seq_logits[i, idx] = 10.0

        # Get gradient
        gradient, current_ll = self.get_gradient(seq_logits, chain=chain)

        # Find positions with largest negative gradient
        # Negative gradient indicates improvement direction
        suggestions = []
        for pos in range(len(sequence)):
            current_aa_idx = seq_indices[pos]
            grad_at_pos = gradient[pos]

            # Find alternative AA with most negative gradient
            for aa_idx, aa in enumerate(self.amino_acids):
                if aa_idx != current_aa_idx:
                    improvement = -grad_at_pos[aa_idx].item()
                    suggestions.append({
                        'position': pos,
                        'from': sequence[pos],
                        'to': aa,
                        'improvement': improvement
                    })

        # Sort by improvement and return top suggestions
        suggestions.sort(key=lambda x: x['improvement'], reverse=True)
        return suggestions[:num_suggestions]