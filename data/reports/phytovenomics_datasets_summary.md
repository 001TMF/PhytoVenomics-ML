# Phytovenomics ML Datasets Summary Report

## Introduction

This report provides a comprehensive overview of all datasets available for the Phytovenomics ML platform, including their structure, content, preprocessing status, and organization for model training.

## Dataset Overview

Total number of datasets: **28**

### Data Distribution

| Dataset Type | Real-World | Synthetic | Total |
| ------------ | ---------- | --------- | ----- |
| Toxins | 589 | 694 | 1283 |
| Antibodies | 277 | 911 | 1188 |
| Binding_Pairs | 10253 | 1974 | 12227 |

## ML-Ready Datasets

These datasets are prepared specifically for machine learning tasks, with proper splits and features.

### Training Datasets

- **antibodies_final_train**
  - Records: 253
  - Features: 23
  - Sources: {'synthetic': 203, 'real': 50}

- **antibodies_train_ml**
  - Records: 298
  - Features: 23
  - Sources: {'synthetic': 239, 'real': 59}

- **binding_final_train**
  - Records: 2969
  - Features: 10
  - Sources: {'real': 2498, 'synthetic': 471}

- **binding_train_ml**
  - Records: 3493
  - Features: 10
  - Sources: {'real': 2927, 'synthetic': 566}

- **toxins_final_train**
  - Records: 255
  - Features: 25
  - Sources: {'synthetic': 128, 'real': 127}

- **toxins_train_ml**
  - Records: 300
  - Features: 25
  - Sources: {'synthetic': 152, 'real': 148}

### Validation Datasets

- **antibodies_validation**
  - Records: 45
  - Features: 23
  - Sources: {'synthetic': 36, 'real': 9}

- **binding_validation**
  - Records: 524
  - Features: 10
  - Sources: {'real': 429, 'synthetic': 95}

- **toxins_validation**
  - Records: 45
  - Features: 25
  - Sources: {'synthetic': 24, 'real': 21}

### Test Datasets

- **antibodies_test_ml**
  - Records: 75
  - Features: 23
  - Sources: {'synthetic': 61, 'real': 14}

- **binding_test_ml**
  - Records: 874
  - Features: 10
  - Sources: {'real': 736, 'synthetic': 138}

- **toxins_test_ml**
  - Records: 76
  - Features: 25
  - Sources: {'synthetic': 48, 'real': 28}

### Cross-Validation Configuration


## Feature Enhancement and Preprocessing

The datasets have undergone various preprocessing steps to enhance their utility for ML models:

### Antibodies_Test

- **Original Dataset Size**: 75 records
- **Added Features**: 15
  - chain_type_code, sheet_fraction, hydrophobicity, instability_index, cysteine_count, isoelectric_point, cdr_score, aromaticity, sequence_complexity, molecular_weight, data_source, turn_fraction, charge_at_pH7, cysteine_percentage, helix_fraction

### Antibodies_Train

- **Original Dataset Size**: 298 records
- **Added Features**: 15
  - chain_type_code, sheet_fraction, hydrophobicity, instability_index, cysteine_count, isoelectric_point, cdr_score, aromaticity, sequence_complexity, molecular_weight, data_source, turn_fraction, charge_at_pH7, cysteine_percentage, helix_fraction

### Toxins_Test

- **Original Dataset Size**: 76 records
- **Added Features**: 16
  - species_code, sheet_fraction, toxin_family_code, hydrophobicity, instability_index, cysteine_count, isoelectric_point, cdr_score, aromaticity, sequence_complexity, molecular_weight, data_source, turn_fraction, charge_at_pH7, cysteine_percentage, helix_fraction

### Toxins_Train

- **Original Dataset Size**: 300 records
- **Added Features**: 16
  - species_code, sheet_fraction, toxin_family_code, hydrophobicity, instability_index, cysteine_count, isoelectric_point, cdr_score, aromaticity, sequence_complexity, molecular_weight, data_source, turn_fraction, charge_at_pH7, cysteine_percentage, helix_fraction


## Dataset File Structure

The datasets are organized in the following directory structure:

```
//
  └── phytovenomics_data_analysis_report.json
  └── phytovenomics_data_report.json
antibody_structures/
  └── sabdab_summary.csv
antibody_structures/human/
  └── antibody_statistics.json
  └── human_antibodies_synthetic.csv
  └── human_antibodies_synthetic.fasta
antibody_structures/real_world/
  └── human_antibodies_metadata.json
  └── human_antibodies_real.csv
  └── human_antibodies_real.fasta
antivenom_research/
  └── antivenom_research.json
ml_ready/
  └── antibodies_balance_stats.json
  └── antibodies_balanced.csv
  └── antibodies_final_train.csv
  └── antibodies_test_ml.csv
  └── antibodies_train_ml.csv
  └── antibodies_validation.csv
  └── antibody_chain_type_stats.json
  └── antibody_cv_folds_chain_type.json
  └── antibody_cv_folds_specificity.json
  └── antibody_specificity_stats.json
  └── binding_binding_type_stats.json
  └── binding_cv_folds_binding_type.json
  └── binding_cv_folds_toxin_type.json
  └── binding_final_train.csv
  └── binding_test_ml.csv
  └── binding_toxin_type_stats.json
  └── binding_train_ml.csv
  └── binding_validation.csv
  └── data_splits_metadata.json
  └── toxin_cv_folds_species.json
  └── toxin_cv_folds_toxin_type.json
  └── toxin_species_stats.json
  └── toxin_toxin_type_stats.json
  └── toxins_balance_stats.json
  └── toxins_balanced.csv
  └── toxins_final_train.csv
  └── toxins_test_ml.csv
  └── toxins_train_ml.csv
  └── toxins_validation.csv
plant_expression/
  └── plant_expression_systems.json
processed/
  └── antibodies_test.csv
  └── antibodies_test.fasta
  └── antibodies_train.csv
  └── antibodies_train.fasta
  └── binding_pairs_test.csv
  └── binding_pairs_train.csv
  └── dataset_metadata.json
  └── toxins_test.csv
  └── toxins_test.fasta
  └── toxins_train.csv
  └── toxins_train.fasta
reports/
  └── phytovenomics_data_report.md
  └── phytovenomics_data_report_summary.json
  └── phytovenomics_enhanced_datasets_report.md
  └── phytovenomics_enhanced_datasets_report_updated.md
snake_toxins/
  └── snake_toxins_processed.csv
  └── snake_toxins_synthetic.csv
  └── snake_toxins_synthetic.fasta
  └── toxin_statistics.json
  └── uniprot_snake_toxins.fasta
  └── uniprot_snake_toxins_metadata.json
snake_toxins/real_world/
  └── snake_toxins_metadata.json
  └── snake_toxins_real.csv
  └── snake_toxins_real.fasta
  └── snake_toxins_structures.json
toxin_antibody_binding/
  └── antivenom_cocktail_suggestion.json
  └── real_binding_pairs.csv
  └── toxin_antibody_binding_pairs.csv
  └── toxin_antibody_binding_statistics.json
```


## Recommendations for ML Model Implementation

- ML datasets have proper training/validation splits ready for model use.
- For imbalanced datasets, consider using the provided class weights or sampling techniques.
- Real-world data is available and should be prioritized for validation.

Additional recommendations:

- **For ESM-2 Transformer Training**: Use the antibody sequence data with CDR annotations for masked language modeling.
- **For RosettaFold Integration**: Leverage the structural data from real-world sources for validation.
- **For Binding Affinity Prediction**: Use the toxin-antibody binding pairs dataset which contains experimental affinity measurements.
- **For Model Validation**: Prioritize real-world data for final validation to ensure biological relevance.


*Report generated on: 2025-05-29 18:24:47*