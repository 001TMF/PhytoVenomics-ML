# Epitope Discovery Engine

## Overview

The Epitope Discovery Engine is a sophisticated component of the Phytovenomics platform designed to identify potential epitopes on snake toxins for targeted antibody design. By leveraging computational prediction algorithms and structural analysis, this module enables the identification of accessible and immunogenic regions on toxin proteins that can be targeted by engineered antibodies.

## Features

- **Sequence-Based Epitope Prediction**: Identifies linear epitopes based on amino acid properties
- **Structure-Based Epitope Mapping**: Determines surface-exposed regions on 3D toxin structures
- **Conserved Epitope Identification**: Finds epitopes conserved across toxin families
- **Epitope Accessibility Analysis**: Evaluates the accessibility of epitopes for antibody binding
- **Immunogenicity Prediction**: Estimates epitope immunogenic potential
- **Cross-Reactivity Assessment**: Evaluates potential cross-reactivity with human proteins

## Architecture

![Epitope Discovery Engine Architecture](../assets/epitope_discovery_architecture.png)

The engine consists of several integrated components:

1. **Sequence Analyzer**: Processes toxin sequences to predict linear epitopes
2. **Structure Analyzer**: Analyzes 3D structures to identify conformational epitopes
3. **Conservation Mapper**: Identifies epitopes conserved across multiple toxins
4. **Epitope Ranker**: Ranks epitopes based on multiple criteria

## Usage

```python
from phytovenomics.epitope_discovery import EpitopeDiscoverer, EpitopeRanker

# Initialize the epitope discoverer
discoverer = EpitopeDiscoverer(config_path="config/epitope_discovery_config.yaml")

# Load toxin sequences
toxin_sequences = ["MKILYLVFIVSLVVGTLEEQGVYGEESRPENTPGSSEPATTEQRGSPGHVCCP...", ...]

# Discover epitopes
epitopes = discoverer.find_epitopes(toxin_sequences)

# Rank epitopes
ranker = EpitopeRanker()
ranked_epitopes = ranker.rank_epitopes(
    epitopes, 
    criteria=["accessibility", "conservation", "immunogenicity"]
)

# Get top epitopes
top_epitopes = ranked_epitopes[:5]

# Export results
discoverer.export_epitopes(top_epitopes, "results/top_epitopes.csv")
```

## Configuration Options

The Epitope Discovery Engine can be configured through the `config/epitope_discovery_config.yaml` file:

```yaml
prediction:
  linear_methods:
    - bepipred
    - abcpred
    - bcepred
  conformational_methods:
    - discotope
    - ellipro

filtering:
  min_length: 8
  max_length: 20
  min_score: 0.7
  exclude_glycosylation_sites: true

ranking:
  weights:
    accessibility: 0.3
    conservation: 0.3
    immunogenicity: 0.2
    surface_exposure: 0.2
```

## Prediction Methods

The engine implements multiple epitope prediction methods:

### Sequence-Based Methods

- **BepiPred**: Uses hidden Markov model for B-cell epitope prediction
- **ABCPred**: Employs artificial neural networks for epitope prediction
- **BCEPred**: Predicts epitopes based on physicochemical properties

### Structure-Based Methods

- **DiscoTope**: Identifies discontinuous epitopes on protein structures
- **ElliPro**: Predicts epitopes based on protrusion index and adjacency
- **PEPOP**: Identifies accessible peptides on protein surfaces

## Epitope Evaluation Criteria

- **Accessibility**: Surface exposure and solvent accessibility
- **Conservation**: Level of conservation across toxin variants
- **Immunogenicity**: Predicted ability to elicit immune response
- **Uniqueness**: Specificity to target toxins (vs. human proteome)
- **Structural Stability**: Likelihood of maintaining native conformation
- **Size and Complexity**: Optimal length and amino acid composition

## Integration

The Epitope Discovery Engine integrates with:

- **Venom Intelligence Module**: Receives toxin sequences and family information
- **Antibody Generator**: Provides ranked epitopes for antibody design
- **Affinity Maturation System**: Provides epitope structures for binding optimization

## API Reference

### EpitopeDiscoverer

```python
class EpitopeDiscoverer:
    def __init__(self, config_path=None):
        """Initialize the epitope discoverer.
        
        Args:
            config_path: Path to configuration file
        """
        pass
        
    def find_epitopes(self, sequences, methods=None):
        """Find epitopes in toxin sequences.
        
        Args:
            sequences: List of toxin sequences
            methods: List of prediction methods to use (default: all configured methods)
            
        Returns:
            DataFrame of predicted epitopes
        """
        pass
        
    def find_conserved_epitopes(self, sequences, min_conservation=0.7):
        """Find epitopes conserved across multiple sequences.
        
        Args:
            sequences: List of toxin sequences
            min_conservation: Minimum conservation score (0-1)
            
        Returns:
            DataFrame of conserved epitopes
        """
        pass
        
    def export_epitopes(self, epitopes, output_path):
        """Export epitope data to file.
        
        Args:
            epitopes: DataFrame of epitopes
            output_path: Path to save results
        """
        pass
```

### EpitopeRanker

```python
class EpitopeRanker:
    def __init__(self, config_path=None):
        """Initialize the epitope ranker.
        
        Args:
            config_path: Path to configuration file
        """
        pass
        
    def rank_epitopes(self, epitopes, criteria=None):
        """Rank epitopes based on specified criteria.
        
        Args:
            epitopes: DataFrame of epitopes
            criteria: List of criteria to use for ranking
            
        Returns:
            DataFrame of ranked epitopes
        """
        pass
        
    def visualize_ranking(self, ranked_epitopes, output_path=None):
        """Visualize epitope ranking.
        
        Args:
            ranked_epitopes: DataFrame of ranked epitopes
            output_path: Optional path to save visualization
        """
        pass
```

## Future Directions

- Integration of deep learning-based epitope prediction
- Molecular dynamics simulation for epitope flexibility analysis
- Expanded cross-reactivity assessment using larger human proteome databases
- Integration of experimental validation feedback

## References

1. Jespersen, M. C., Peters, B., Nielsen, M., & Marcatili, P. (2017). BepiPred-2.0: improving sequence-based B-cell epitope prediction using conformational epitopes. Nucleic acids research, 45(W1), W24-W29.

2. Kringelum, J. V., Lundegaard, C., Lund, O., & Nielsen, M. (2012). Reliable B cell epitope predictions: impacts of method development and improved benchmarking. PLoS computational biology, 8(12), e1002829.

3. Saha, S., & Raghava, G. P. S. (2006). Prediction of continuous B‐cell epitopes in an antigen using recurrent neural network. Proteins: Structure, Function, and Bioinformatics, 65(1), 40-48.

4. Kringelum, J. V., Nielsen, M., Padkjær, S. B., & Lund, O. (2013). Structural analysis of B-cell epitopes in antibody:protein complexes. Molecular immunology, 53(1-2), 24-34.