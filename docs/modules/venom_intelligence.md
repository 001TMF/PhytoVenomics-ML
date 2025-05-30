# Venom Intelligence Module

## Overview

The Venom Intelligence Module is a core component of the Phytovenomics platform that analyzes and clusters snake toxins using a comprehensive database. This module processes toxin sequences, identifies toxin families, and provides insights into venom composition across different snake species.

## Features

- **Toxin Classification**: Categorizes toxins into families (e.g., phospholipases, metalloproteases, three-finger toxins)
- **Sequence Analysis**: Performs multiple sequence alignment and evolutionary analysis
- **Clustering**: Groups toxins based on sequence similarity and functional properties
- **Cross-Species Analysis**: Compares venom compositions across different snake species
- **Toxin Prioritization**: Identifies high-priority toxins for antibody development

## Architecture

![Venom Intelligence Module Architecture](../assets/venom_intelligence_architecture.png)

The module consists of several components:

1. **Toxin Database Manager**: Manages the database of toxin sequences and metadata
2. **Sequence Analyzer**: Performs sequence analysis and alignment
3. **Toxin Clusterer**: Groups toxins based on various criteria
4. **Evolutionary Analyzer**: Analyzes evolutionary relationships between toxins

## Usage

```python
from phytovenomics.venom_data import ToxinDatabase
from phytovenomics.venom_intelligence import ToxinClusterer, EvolutionaryAnalyzer

# Load toxin database
toxin_db = ToxinDatabase()
toxin_db.load_data("data/snake_toxins/snake_toxins_processed.csv")

# Cluster toxins
clusterer = ToxinClusterer()
clusters = clusterer.cluster_toxins(toxin_db.get_sequences(), method="hierarchical", n_clusters=10)

# Analyze evolutionary relationships
evol_analyzer = EvolutionaryAnalyzer()
phylo_tree = evol_analyzer.build_phylogenetic_tree(toxin_db.get_sequences())
evol_analyzer.visualize_tree(phylo_tree, "results/toxin_phylogeny.png")

# Get high-priority toxins
priority_toxins = toxin_db.get_priority_toxins(species="Bothrops asper", top_n=5)
```

## Configuration Options

The Venom Intelligence Module can be configured through the `config/venom_intelligence_config.yaml` file:

```yaml
database:
  update_frequency: monthly
  sources:
    - UniProt
    - VenomZone
    - ToxProt

clustering:
  default_method: hierarchical
  linkage: average
  distance_metric: blosum62

phylogeny:
  bootstrap_iterations: 100
  substitution_model: JTT
```

## Data Sources

The module integrates data from multiple sources:

- **UniProt**: Comprehensive protein database including toxin sequences
- **VenomZone**: Specialized database for venom proteins
- **ToxProt**: Annotated toxin protein database
- **Published literature**: Curated toxin sequences from research papers

## Integration

The Venom Intelligence Module integrates with:

- **Epitope Discovery Engine**: Provides toxin sequences for epitope identification
- **Antibody Generator**: Supplies priority toxins for antibody design
- **Cocktail Formulator**: Informs cocktail composition based on toxin clustering

## API Reference

### ToxinDatabase

```python
class ToxinDatabase:
    def __init__(self, config_path=None):
        """Initialize the toxin database.
        
        Args:
            config_path: Path to configuration file
        """
        pass
        
    def load_data(self, data_path):
        """Load toxin data from file.
        
        Args:
            data_path: Path to toxin data file
        """
        pass
        
    def get_sequences(self, filters=None):
        """Get toxin sequences with optional filtering.
        
        Args:
            filters: Dictionary of filters (e.g., species, family)
            
        Returns:
            DataFrame of toxin sequences
        """
        pass
        
    def get_priority_toxins(self, species=None, top_n=10):
        """Get high-priority toxins for antibody development.
        
        Args:
            species: Optional species filter
            top_n: Number of toxins to return
            
        Returns:
            List of high-priority toxin sequences
        """
        pass
```

### ToxinClusterer

```python
class ToxinClusterer:
    def __init__(self, config_path=None):
        """Initialize the toxin clusterer.
        
        Args:
            config_path: Path to configuration file
        """
        pass
        
    def cluster_toxins(self, sequences, method="hierarchical", n_clusters=None):
        """Cluster toxins based on sequence similarity.
        
        Args:
            sequences: DataFrame of toxin sequences
            method: Clustering method ("hierarchical", "kmeans", "dbscan")
            n_clusters: Number of clusters (for applicable methods)
            
        Returns:
            DataFrame with cluster assignments
        """
        pass
        
    def visualize_clusters(self, clusters, output_path=None):
        """Visualize toxin clusters.
        
        Args:
            clusters: Clustering results
            output_path: Optional path to save visualization
        """
        pass
```

## Future Directions

- Integration of 3D structure prediction for toxins
- Machine learning-based toxicity prediction
- Expanded cross-reactivity analysis
- Integration of mass spectrometry data for improved toxin identification

## References

1. Lomonte, B., & Calvete, J. J. (2017). Strategies in 'snake venomics' aiming at an integrative view of compositional, functional, and immunological characteristics of venoms. Journal of Venomous Animals and Toxins including Tropical Diseases, 23.

2. Tasoulis, T., & Isbister, G. K. (2017). A review and database of snake venom proteomes. Toxins, 9(9), 290.

3. Ainsworth, S., et al. (2018). The paraspecific neutralisation of snake venom induced coagulopathy by antivenoms. Communications biology, 1(1), 1-14.