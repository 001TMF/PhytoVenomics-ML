# Cocktail Formulator

## Overview

The Cocktail Formulator is a specialized component of the Phytovenomics platform responsible for designing optimal antibody cocktails that effectively neutralize snake venoms. This module leverages computational methods to combine individual antibodies into synergistic mixtures, maximizing coverage against diverse toxins while minimizing the total number of components required.

## Features

- **Coverage Optimization**: Maximizes toxin neutralization coverage across snake species
- **Synergy Analysis**: Identifies antibody combinations with synergistic effects
- **Component Minimization**: Minimizes the number of antibodies needed for effective neutralization
- **Cross-Reactivity Assessment**: Evaluates antibody cross-reactivity across toxin families
- **Toxin Prioritization**: Weights cocktail composition based on toxin medical importance
- **Geographic Tailoring**: Designs region-specific cocktails based on local snake species
- **Manufacturing Constraints Integration**: Incorporates expression yield and production considerations

## Architecture

![Cocktail Formulator Architecture](../assets/cocktail_formulator_architecture.png)

The Cocktail Formulator consists of several integrated components:

1. **Coverage Analyzer**: Analyzes toxin coverage of individual antibodies
2. **Synergy Calculator**: Assesses synergistic effects between antibody pairs
3. **Cocktail Optimizer**: Optimizes antibody combinations using mathematical programming
4. **Production Evaluator**: Assesses manufacturability of proposed cocktails
5. **Stability Predictor**: Evaluates the stability of antibody combinations

## Usage

```python
from phytovenomics.cocktail_formulation import CocktailFormulator

# Initialize the cocktail formulator
formulator = CocktailFormulator(config_path="config/cocktail_formulator_config.yaml")

# Load antibodies and their toxin neutralization profiles
antibodies = [...] # List of antibody objects with neutralization data
toxin_map = {...} # Dictionary mapping toxins to their medical importance scores

# Create optimized cocktail formulations
cocktails = formulator.design_cocktails(
    antibodies=antibodies,
    toxin_importance=toxin_map,
    max_components=8,
    coverage_threshold=0.95
)

# Evaluate cocktail performance
performance = formulator.evaluate_cocktails(cocktails, detailed=True)

# Export results
formulator.export_cocktails(cocktails, performance, "results/optimized_cocktails.csv")

# Visualize cocktail coverage
formulator.visualize_coverage(cocktails[0], "results/cocktail_coverage.png")
```

## Configuration Options

The Cocktail Formulator can be configured through the `config/cocktail_formulator_config.yaml` file:

```yaml
optimization:
  algorithm: "mixed_integer_programming"
  objective: "minimize_components"
  coverage_constraint: 0.9
  max_components: 10
  time_limit_seconds: 600

toxin_weighting:
  method: "medical_importance"
  default_weight: 1.0
  critical_factor: 5.0

geographic_tailoring:
  enabled: true
  regions:
    - name: "Sub-Saharan Africa"
      species:
        - "Bitis arietans"
        - "Dendroaspis polylepis"
        - "Naja nigricollis"
    - name: "South Asia"
      species:
        - "Naja naja"
        - "Daboia russelii"
        - "Bungarus caeruleus"
        
manufacturing_constraints:
  min_expression_threshold: 0.5  # g/L
  compatibility_matrix_path: "data/antibody_compatibility_matrix.csv"
  cost_model: "linear"
```

## Optimization Methodology

The Cocktail Formulator employs several optimization approaches:

### Mathematical Programming

The core optimization uses mixed-integer linear programming (MILP) to find optimal antibody combinations:

- **Objective Function**: Minimize the weighted sum of included antibodies
- **Constraints**:
  - Achieve required coverage for each toxin
  - Respect maximum component count
  - Consider antibody compatibility constraints
- **Decision Variables**: Binary variables indicating antibody inclusion

### Greedy Approximation

For larger problems, a greedy algorithm approach is also available:

1. Start with an empty cocktail
2. Iteratively add the antibody providing the greatest increase in coverage
3. Continue until coverage threshold or component limit is reached

### Genetic Algorithm

For complex cases with non-linear synergy effects:

1. Encode cocktail compositions as binary strings
2. Evaluate fitness based on coverage and component count
3. Evolve population through selection, crossover, and mutation
4. Extract best solutions after convergence

## Coverage Calculation

The system employs several methods to calculate toxin coverage:

- **Binary Coverage**: Simple binding/non-binding determination
- **Affinity-Based**: Coverage proportional to binding affinity
- **Neutralization-Based**: Based on in silico neutralization predictions
- **Synergy-Aware**: Includes cooperative effects between antibodies

## Geographic Tailoring

The Cocktail Formulator supports geographic customization:

- **Region Definition**: Define regions based on snake species presence
- **Species Importance**: Weight species by medical importance in each region
- **Venom Composition**: Account for regional variation in venom composition
- **Custom Formulations**: Generate specialized formulations for different regions

## Integration

The Cocktail Formulator integrates with:

- **Venom Intelligence Module**: Receives toxin prioritization data
- **Antibody Generator**: Receives candidate antibodies
- **Affinity Maturation System**: Receives binding affinity data
- **Production Optimization System**: Feeds cocktail formulations for production planning

## API Reference

### CocktailFormulator

```python
class CocktailFormulator:
    def __init__(self, config_path=None):
        """Initialize the cocktail formulator.
        
        Args:
            config_path: Path to configuration file
        """
        pass
        
    def design_cocktails(self, antibodies, toxin_importance=None, max_components=10, coverage_threshold=0.9):
        """Design optimized antibody cocktails.
        
        Args:
            antibodies: List of antibody objects with neutralization data
            toxin_importance: Dictionary mapping toxins to importance scores
            max_components: Maximum number of antibodies in a cocktail
            coverage_threshold: Minimum required toxin coverage
            
        Returns:
            List of cocktail formulations
        """
        pass
        
    def evaluate_cocktail(self, cocktail, detailed=False):
        """Evaluate a cocktail's performance.
        
        Args:
            cocktail: Cocktail formulation
            detailed: Whether to return detailed metrics
            
        Returns:
            Dictionary of performance metrics
        """
        pass
        
    def evaluate_cocktails(self, cocktails, detailed=False):
        """Evaluate multiple cocktails' performance.
        
        Args:
            cocktails: List of cocktail formulations
            detailed: Whether to return detailed metrics
            
        Returns:
            List of performance metric dictionaries
        """
        pass
        
    def calculate_synergy(self, antibody1, antibody2):
        """Calculate synergy between two antibodies.
        
        Args:
            antibody1: First antibody
            antibody2: Second antibody
            
        Returns:
            Synergy score
        """
        pass
        
    def create_regional_cocktail(self, antibodies, region_name):
        """Create a cocktail tailored to a specific geographic region.
        
        Args:
            antibodies: List of available antibodies
            region_name: Name of the region
            
        Returns:
            Regionally optimized cocktail
        """
        pass
        
    def export_cocktails(self, cocktails, performance_data, output_path):
        """Export cocktail formulations to file.
        
        Args:
            cocktails: List of cocktail formulations
            performance_data: Performance metrics for each cocktail
            output_path: Path to save results
        """
        pass
        
    def visualize_coverage(self, cocktail, output_path=None):
        """Visualize a cocktail's toxin coverage.
        
        Args:
            cocktail: Cocktail formulation
            output_path: Optional path to save visualization
        """
        pass
```

### CocktailOptimizer

```python
class CocktailOptimizer:
    def __init__(self, algorithm="milp", config_path=None):
        """Initialize the cocktail optimizer.
        
        Args:
            algorithm: Optimization algorithm ("milp", "greedy", "genetic")
            config_path: Path to configuration file
        """
        pass
        
    def optimize(self, coverage_matrix, toxin_weights=None, max_components=10, coverage_threshold=0.9):
        """Run optimization to find optimal antibody combinations.
        
        Args:
            coverage_matrix: Matrix of antibody-toxin coverage values
            toxin_weights: Optional weights for toxins
            max_components: Maximum number of antibodies in the cocktail
            coverage_threshold: Minimum required toxin coverage
            
        Returns:
            Optimal antibody combination
        """
        pass
```

## Future Directions

- Integration of pharmacokinetic modeling for sustained protection
- Machine learning prediction of antibody-antibody interactions
- Multi-objective optimization considering cost, shelf stability, and coverage
- Integration of experimental neutralization data for improved accuracy
- Development of modular cocktail systems with mix-and-match components

## References

1. Laustsen, A. H., Engmark, M., Milbo, C., Johannesen, J., Lomonte, B., Gutiérrez, J. M., & Lohse, B. (2016). From fangs to pharmacology: The future of snakebite envenoming therapy. Current Pharmaceutical Design, 22(34), 5270-5293.

2. Ratanabanangkoon, K., et al. (2018). A novel in vitro potency assay of antisera against Thai Naja kaouthia based on nicotinic acetylcholine receptor binding. Scientific Reports, 8(1), 1-9.

3. Ainsworth, S., et al. (2018). The paraspecific neutralisation of snake venom induced coagulopathy by antivenoms. Communications Biology, 1(1), 1-14.

4. Jenkins, T. P., et al. (2019). Toxivenomics and antivenom profiling of the Eastern green mamba snake (Dendroaspis angusticeps). Journal of Proteomics, 181, 183-192.

5. Laustsen, A. H. (2018). Toxin-centric development approach for next-generation antivenoms. Toxicon, 150, 195-197.