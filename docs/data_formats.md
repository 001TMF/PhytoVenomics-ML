# Phytovenomics Data Formats

## Overview

The Phytovenomics platform uses standardized data formats for consistent data handling across modules. This document describes the key data structures, file formats, and conventions used throughout the platform to ensure interoperability and efficient data flow between components.

## Data Directory Structure

The platform organizes data in a hierarchical directory structure:

```
data/
├── snake_toxins/           # Snake toxin sequence and structure data
├── antibodies/             # Generated and reference antibody data
├── epitopes/               # Identified epitopes and binding sites
├── cocktails/              # Antibody cocktail formulations
├── structures/             # 3D structure files
├── validation/             # Validation results and metrics
├── production/             # Production parameters and yields
└── reference/              # Reference databases and lookups
```

## Key Data Formats

### Toxin Database Format

Toxin data is stored in CSV format with the following fields:

```csv
toxin_id,species,family,subfamily,sequence,molecular_weight,isoelectric_point,lethality_score,medical_importance,structure_available
```

Example:
```
toxin_id,species,family,subfamily,sequence,molecular_weight,isoelectric_point,lethality_score,medical_importance,structure_available
TX001,Naja naja,3FTX,Short neurotoxin,MKTLLLTLVVVTIVCLDLGYTLECHNQQSSQPPTTKTCSGETNCYKKWWSDHRGTIIERGCGCPTVKPGIKLNCCTTDKCNN,7825.8,8.3,4.7,5,TRUE
TX002,Bothrops jararaca,SVMP,P-I,MIQVLLVTICLAVFPYQGSSIILESGNVNDYEVVYPQKVTALPKGAVQPKYEDAMQYEFKVNGESP...,22635.4,5.2,3.2,4,FALSE
```

### Antibody Sequence Format

Antibody sequences are stored in FASTA format with structured headers:

```
>AB001|IgG|Heavy|Phytovenomics|Target:TX001|Affinity:0.5nM|Plant:Optimized
EVQLVESGGGLVQPGGSLRLSCAASSFTFSSYWMSWVRQAPGKGLEWVANIKQDGSEKYYVDSVKGRFTISRDNAKNSLYLQMNSLRDEDTAVYYCAREGHTFDIWGQGTMVTVSS
>AB001|IgG|Light|Phytovenomics|Target:TX001|Affinity:0.5nM|Plant:Optimized
EIVLTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYDASSRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQYGSSPLTFGGGTKVEIK
```

### Epitope Data Format

Epitope data is stored in JSON format:

```json
{
  "epitope_id": "EP001",
  "toxin_id": "TX001",
  "start_position": 25,
  "end_position": 38,
  "sequence": "NQQSSQPPTTKTCS",
  "score": 0.87,
  "accessibility": 0.92,
  "conservation": 0.78,
  "structure_coordinates": [
    {"residue": "N", "position": 25, "x": 10.4, "y": 5.3, "z": 1.2},
    {"residue": "Q", "position": 26, "x": 9.8, "y": 4.9, "z": 1.8},
    "..."  
  ],
  "properties": {
    "hydrophobicity": -0.5,
    "charge": 1,
    "flexibility": 0.6
  },
  "targeting_antibodies": ["AB001", "AB004", "AB009"]
}
```

### Structure Data Format

Protein structures are stored in standard PDB format for atomic coordinates and mmCIF format for complex structures. Simplified example of PDB format:

```
ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00 37.55      A    N
ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 37.80      A    C
ATOM      3  C   MET A   1      26.913  26.639   3.531  1.00 37.65      A    C
ATOM      4  O   MET A   1      27.886  26.463   4.263  1.00 38.18      A    O
ATOM      5  CB  MET A   1      25.112  24.880   3.649  1.00 38.31      A    C
...
```

### Binding Prediction Data Format

Binding prediction results are stored in JSON format:

```json
{
  "antibody_id": "AB001",
  "toxin_id": "TX001",
  "binding_energy": -12.3,
  "confidence": 0.88,
  "interface_area": 1240.5,
  "hydrogen_bonds": 12,
  "salt_bridges": 3,
  "interaction_residues": {
    "antibody": ["H:Y32", "H:W52", "H:Y56", "H:R99", "L:Y32", "L:N92", "L:W96"],
    "toxin": ["N25", "Q26", "S28", "P30", "T33", "K34", "T35"]
  },
  "prediction_method": "hybrid",
  "timestamp": "2025-05-01T14:22:33Z"
}
```

### Cocktail Formulation Format

Antibody cocktail formulations are stored in JSON format:

```json
{
  "cocktail_id": "COC001",
  "name": "PanAfrican-1",
  "description": "Broad-spectrum cocktail optimized for African elapid venoms",
  "components": [
    {"antibody_id": "AB001", "concentration": 5.0, "units": "mg/ml"},
    {"antibody_id": "AB007", "concentration": 3.5, "units": "mg/ml"},
    {"antibody_id": "AB012", "concentration": 4.0, "units": "mg/ml"}
  ],
  "target_species": ["Naja melanoleuca", "Dendroaspis polylepis", "Naja nigricollis"],
  "coverage_metrics": {
    "overall": 0.92,
    "by_species": {
      "Naja melanoleuca": 0.95,
      "Dendroaspis polylepis": 0.88,
      "Naja nigricollis": 0.94
    },
    "by_toxin_family": {
      "3FTX": 0.97,
      "PLA2": 0.85,
      "SVMP": 0.65
    }
  },
  "validation_status": "passed",
  "production_parameters": {
    "expression_system": "Nicotiana benthamiana",
    "expected_yield": 120,
    "units": "mg/kg fresh weight"
  }
}
```

### Validation Metrics Format

Validation metrics are stored in JSON format:

```json
{
  "validation_id": "VAL001",
  "entity_type": "antibody",
  "entity_id": "AB001",
  "timestamp": "2025-05-02T09:14:22Z",
  "validation_level": "comprehensive",
  "overall_status": "passed",
  "overall_score": 0.89,
  "metrics": {
    "sequence_quality": {
      "score": 0.95,
      "humanness": 0.91,
      "unusual_residues": 0,
      "status": "passed"
    },
    "developability": {
      "score": 0.88,
      "hydrophobicity_index": 0.7,
      "charge_at_ph7": 1.2,
      "aggregation_score": 18,
      "status": "passed"
    },
    "predicted_binding": {
      "score": 0.92,
      "binding_energy": -11.4,
      "interface_area": 825,
      "shape_complementarity": 0.78,
      "status": "passed"
    },
    "predicted_structure": {
      "score": 0.85,
      "packing_quality": 0.83,
      "ramachandran_outliers": 1.2,
      "status": "passed"
    },
    "cross_reactivity": {
      "score": 0.98,
      "human_proteome_hits": 0,
      "off_target_binding": 0.12,
      "status": "passed"
    },
    "plant_expressibility": {
      "score": 0.79,
      "predicted_yield": 85,
      "proteolysis_sites": 2,
      "status": "passed"
    }
  },
  "issues": [],
  "recommendations": []
}
```

### Experimental Data Format

Experimental validation data is stored in CSV format:

```csv
antibody_id,toxin_id,assay_type,assay_date,parameter,value,units,operator,platform,notes
AB001,TX001,binding_kinetics,2025-04-15,kon,3.2e5,1/Ms,,Octet RED96,
AB001,TX001,binding_kinetics,2025-04-15,koff,1.6e-4,1/s,,Octet RED96,
AB001,TX001,binding_kinetics,2025-04-15,KD,0.5,nM,,Octet RED96,calculated
AB001,TX001,neutralization,2025-04-20,IC50,2.3,nM,,cell-based,HEK293 cells
AB001,TX001,thermal_stability,2025-04-18,Tm,72.5,C,,DSF,PBS pH 7.4
```

### Production Data Format

Production data is stored in JSON format:

```json
{
  "production_id": "PROD001",
  "entity_id": "AB001",
  "entity_type": "antibody",
  "expression_system": "Nicotiana benthamiana",
  "construct_id": "CONS001",
  "batch": "B2025-05-01",
  "scale": "pilot",
  "parameters": {
    "infiltration_OD": 0.6,
    "growth_temperature": 24,
    "light_cycle": [16, 8],
    "harvest_day": 6,
    "extraction_buffer": "PBS pH 7.4 + 10mM ascorbic acid"
  },
  "results": {
    "biomass": 5.2,
    "biomass_units": "kg",
    "total_yield": 468,
    "yield_units": "mg",
    "specific_yield": 90,
    "specific_yield_units": "mg/kg FW",
    "purity": 98.5,
    "purity_units": "%",
    "activity_retention": 94,
    "activity_retention_units": "%"
  },
  "quality": {
    "aggregation": 2.1,
    "aggregation_units": "%",
    "endotoxin": 0.05,
    "endotoxin_units": "EU/mg",
    "host_cell_protein": 12,
    "host_cell_protein_units": "ppm"
  },
  "timestamp": "2025-05-10T16:45:33Z",
  "operator": "J. Smith"
}
```

## Data Interchange Formats

### Module Input/Output Format

Modules exchange data using a standardized JSON format:

```json
{
  "request_id": "REQ123456",
  "module": "epitope_discovery",
  "action": "find_epitopes",
  "parameters": {
    "toxin_ids": ["TX001", "TX002", "TX003"],
    "min_score": 0.7,
    "conservation_required": true
  },
  "data": {
    "toxins": [
      {"toxin_id": "TX001", "sequence": "MKTLLLTLVVVTIVCLDLGYT..."},
      {"toxin_id": "TX002", "sequence": "MIQVLLVTICLAVFPYQGSSI..."},
      {"toxin_id": "TX003", "sequence": "MIAFIVMIVVLMTAWEYPSYG..."}
    ]
  },
  "timestamp": "2025-05-15T10:23:45Z"
}
```

### Analysis Results Format

Standardized JSON format for analysis results:

```json
{
  "response_id": "RES123456",
  "request_id": "REQ123456",
  "module": "epitope_discovery",
  "action": "find_epitopes",
  "status": "success",
  "results": {
    "epitopes": [
      {
        "epitope_id": "EP001",
        "toxin_id": "TX001",
        "start_position": 25,
        "end_position": 38,
        "sequence": "NQQSSQPPTTKTCS",
        "score": 0.87
      },
      {
        "epitope_id": "EP002",
        "toxin_id": "TX002",
        "start_position": 42,
        "end_position": 53,
        "sequence": "KVTALPKGAVQP",
        "score": 0.81
      }
    ],
    "metadata": {
      "algorithm_version": "2.3.0",
      "parameters_used": {
        "min_score": 0.7,
        "conservation_required": true
      },
      "execution_time": 12.4
    }
  },
  "timestamp": "2025-05-15T10:24:12Z"
}
```

## File Naming Conventions

The platform follows consistent file naming conventions:

- **Snake toxin files**: `toxin_[database]_[date].[format]` (e.g., `toxin_uniprot_20250301.csv`)
- **Antibody files**: `antibody_[id]_[chain].[format]` (e.g., `antibody_AB001_heavy.fasta`)
- **Epitope files**: `epitope_[toxin_id]_[date].[format]` (e.g., `epitope_TX001_20250415.json`)
- **Structure files**: `structure_[entity_id]_[type].[format]` (e.g., `structure_AB001_complex.pdb`)
- **Validation files**: `validation_[entity_id]_[date].[format]` (e.g., `validation_AB001_20250420.json`)
- **Cocktail files**: `cocktail_[id]_[version].[format]` (e.g., `cocktail_COC001_v2.json`)
- **Production files**: `production_[entity_id]_[batch].[format]` (e.g., `production_AB001_B20250501.json`)

## Version Control

Data versioning follows these conventions:

- **Major database updates**: Increment first number (e.g., v2.0.0 to v3.0.0)
- **Addition of new data fields**: Increment second number (e.g., v2.0.0 to v2.1.0)
- **Data corrections/updates**: Increment third number (e.g., v2.1.0 to v2.1.1)
- **Date stamps**: Included in filenames for temporal tracking (YYYYMMDD format)

## Data Flow Between Modules

The following diagram illustrates the data flow between major platform modules:

```
┌──────────────────┐     toxin_data.csv     ┌───────────────────┐
│                  │─────────────────────────▶                   │
│ Venom            │                         │ Epitope           │
│ Intelligence     │                         │ Discovery         │
│ Module           │◀────────────────────────│ Engine            │
│                  │     epitopes.json       │                   │
└──────────────────┘                         └───────────────────┘
         │                                             │
         │                                             │
         │                                             │
         ▼                                             ▼
┌──────────────────┐    antibody_design.json  ┌───────────────────┐
│                  │◀────────────────────────│                   │
│ Hybrid           │                         │ Antibody          │
│ Prediction       │                         │ Generator         │
│ Pipeline         │─────────────────────────▶                   │
│                  │  binding_predictions.json│                   │
└──────────────────┘                         └───────────────────┘
         │                                             │
         │                                             │
         │                                             │
         ▼                                             ▼
┌──────────────────┐    optimized_antibody.json┌───────────────────┐
│                  │◀────────────────────────│                   │
│ Validation       │                         │ Affinity          │
│ Metrics          │                         │ Maturation        │
│ Tracker          │─────────────────────────▶                   │
│                  │   validation_metrics.json│                   │
└──────────────────┘                         └───────────────────┘
         │                                             │
         │                                             │
         │                                             │
         ▼                                             ▼
┌──────────────────┐     cocktail_design.json  ┌───────────────────┐
│                  │◀────────────────────────│                   │
│ Plant            │                         │ Cocktail          │
│ Production       │                         │ Formulator        │
│ Optimization     │─────────────────────────▶                   │
│                  │   production_metrics.json│                   │
└──────────────────┘                         └───────────────────┘
```

## Data Quality Control

The platform implements several data quality control measures:

1. **Schema Validation**: All data files are validated against defined schemas
2. **Range Checks**: Numeric values are verified to be within expected ranges
3. **Consistency Checks**: Cross-references between data files are verified
4. **Completeness Checks**: Required fields are verified to be present and populated
5. **Format Checks**: Data formats are verified to match expectations

## Data Security and Privacy

The platform implements the following data security measures:

1. **Access Control**: Role-based access control for data files
2. **Encryption**: Sensitive data is encrypted at rest
3. **Audit Logs**: All data access and modifications are logged
4. **Backup**: Regular automated backups of all data
5. **Anonymization**: Personal information is anonymized where applicable

## Future Data Format Extensions

The platform is designed to accommodate future extensions:

1. **Clinical Data Integration**: Formats for clinical trial and patient data
2. **Manufacturing Scale-up**: Extended production data for commercial scale
3. **Regulatory Submissions**: Formats aligned with regulatory requirements
4. **Multi-omics Integration**: Integration with genomics, proteomics data
5. **Field Testing Data**: Mobile data collection for field validation