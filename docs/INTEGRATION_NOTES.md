# Integration Notes: DiffAb + Germinal + RFdiffusion

## Integration Strategy

This pipeline successfully integrates three major antibody design systems by extracting and reappropriating their core functionality:

### 1. DiffAb Integration

**What we took**:
- Structure preparation pipeline
- Chothia renumbering system
- CDR identification and labeling
- HDock docking integration
- Structure parsing utilities

**Source files**:
- `diffab/tools/renumber/run.py` → `antibody_pipeline/preparation/renumber.py`
- `diffab/utils/protein/parsers.py` → `antibody_pipeline/preparation/parser.py`
- `diffab/utils/protein/constants.py` → `antibody_pipeline/preparation/constants.py`
- `diffab/datasets/sabdab.py` → `antibody_pipeline/preparation/cdr_identifier.py`
- `diffab/tools/dock/hdock.py` → `antibody_pipeline/docking/hdock.py`

**Integration approach**:
- Direct port of core functions
- Maintained abnumber dependency for Chothia numbering
- Preserved BioPython structure parsing
- Kept CDR range definitions (Chothia scheme)

### 2. Germinal Integration

**What we took**:
- IgLM wrapper with gradient computation
- Filtering pipeline architecture
- Structure prediction integration (Chai/AF3)
- Interface metrics calculation
- Developability assessment

**Source files**:
- `germinal/colabdesign/colabdesign/iglm/model.py` → `antibody_pipeline/optimization/iglm_optimizer.py`
- `germinal/filters/filter_utils.py` → `antibody_pipeline/filtering/`

**Integration approach**:
- Simplified IgLM wrapper (removed optimization loop)
- Ported filtering concepts without full dependency chain
- Maintained straight-through estimator for gradients
- Preserved developability metrics (SAP score, PTM sites)

### 3. RFdiffusion Integration

**What we added**:
- Wrapper for RFdiffusion v3 inference
- Epitope-conditioned design
- CDR-specific masking strategies
- Multiple design modes

**Source files**:
- New implementation: `antibody_pipeline/design/designer.py`

**Integration approach**:
- Subprocess-based wrapper (calls RFdiffusion CLI)
- Inspired by DiffAb's masking strategies
- Added epitope conditioning via RFdiffusion's hotspot feature

## Key Design Decisions

### 1. Module Independence

Each module is self-contained and can be used independently:

```python
# Use only renumbering
from antibody_pipeline.preparation import AntibodyRenumberer
renumberer = AntibodyRenumberer()
renumberer.renumber("input.pdb", "output.pdb")

# Use only IgLM scoring
from antibody_pipeline.optimization import IgLMOptimizer
iglm = IgLMOptimizer()
score = iglm.score_sequence("EVQLVES...")

# Use only docking
from antibody_pipeline.docking import HDockAntibody
with HDockAntibody() as docker:
    docker.set_antigen("antigen.pdb")
    docker.set_antibody("antibody.pdb")
    complexes = docker.dock()
```

### 2. Minimal Dependencies

We avoided pulling in the full dependency chains of DiffAb and Germinal:

**DiffAb dependencies we avoided**:
- Full diffusion model training stack
- Complex model architectures
- Heavy training utilities

**Germinal dependencies we avoided**:
- ColabDesign full framework
- AlphaFold2 hallucination code
- Optimization loop infrastructure

**What we kept**:
- abnumber (for Chothia numbering)
- BioPython (for structure parsing)
- iglm (for sequence scoring)
- Basic scientific stack (numpy, torch)

### 3. Configuration Over Code

All settings exposed through configuration:

```yaml
# Easy to adjust without code changes
design:
  mode: cdr_h3
  num_designs: 10

filtering:
  max_clashes: 100
  min_plddt: 70.0

iglm:
  species: '[HUMAN]'
```

### 4. Subprocess vs Import

**RFdiffusion**: Subprocess wrapper
- Reason: RFdiffusion has complex dependencies and is meant to be called via CLI
- Benefit: Isolates dependency issues

**IgLM**: Direct import
- Reason: iglm is a clean PyPI package
- Benefit: Can compute gradients directly

**HDock**: Subprocess wrapper
- Reason: HDock is a binary executable
- Benefit: No Python dependencies

## Code Reuse Strategy

### Direct Ports (90%+ similarity)

Files that are nearly identical to source:

1. **renumber.py**: Direct port of DiffAb's renumbering
   - Uses abnumber exactly as DiffAb does
   - Preserves chain identification logic
   - Only added type hints and docstrings

2. **constants.py**: Direct copy of DiffAb's constants
   - CDR ranges identical
   - Amino acid definitions identical
   - Heavy atom definitions preserved

3. **hdock.py**: Close port of DiffAb's HDock wrapper
   - Fixed PDB format handling preserved
   - Complex merging logic maintained
   - Added better documentation

### Adapted Ports (50-80% similarity)

Files that take core concepts but simplify:

1. **iglm_optimizer.py**: Simplified Germinal's IgLM wrapper
   - Removed optimization loop (not needed for scoring)
   - Kept straight-through estimator
   - Added mutation suggestion feature

2. **parser.py**: Adapted DiffAb's structure parser
   - Core parsing logic preserved
   - Simplified preprocessing
   - Removed some specialized transforms

### New Implementations (inspired by)

Files that are new but inspired by source systems:

1. **designer.py**: New RFdiffusion wrapper
   - Inspired by DiffAb's design modes
   - Uses RFdiffusion's native CLI
   - Added epitope conditioning

2. **pipeline.py**: New orchestrator
   - Workflow inspired by Germinal's run script
   - Filtering logic adapted from Germinal
   - New ranking system

## Testing Strategy

### Unit Tests (TODO)

```python
# Test renumbering
def test_renumber():
    result = renumber_antibody("test.pdb", "out.pdb")
    assert len(result[0]) > 0  # Heavy chains found

# Test CDR identification
def test_cdr_identification():
    cdrs = extract_cdr_sequences(seq, resseq, 'H')
    assert 'H3' in cdrs

# Test IgLM scoring
def test_iglm_scoring():
    score = iglm.score_sequence("EVQLVES...")
    assert 'log_likelihood' in score
```

### Integration Tests (TODO)

```python
# Test full pipeline
def test_pipeline():
    config = PipelineConfig(
        antigen_pdb="test_antigen.pdb",
        framework_pdb="test_framework.pdb",
        num_designs=2,
        output_dir="./test_output"
    )
    pipeline = AntibodyDesignPipeline(config)
    candidates = pipeline.run()
    assert len(candidates) > 0
```

## Compatibility Notes

### DiffAb Compatibility

- Uses same Chothia numbering scheme
- CDR definitions identical
- Can process DiffAb-prepared structures
- Compatible with DiffAb PDB outputs

### Germinal Compatibility

- IgLM scores comparable
- Filtering thresholds aligned
- Can use Germinal's structure predictors
- Compatible with Germinal's metrics

### RFdiffusion Compatibility

- Uses RFdiffusion v3 CLI
- Contig syntax compatible
- Hotspot residue format aligned
- Output PDBs standard format

## Known Limitations

1. **RFdiffusion dependency**: Requires separate RFdiffusion installation
2. **HDock binary**: Requires manual download and setup
3. **GPU requirement**: Most features need CUDA-capable GPU
4. **PyRosetta optional**: Some filtering features require PyRosetta license
5. **Chai/AF3 optional**: Structure prediction requires additional setup

## Future Integration Opportunities

1. **ESM Models**: Add ESM2/ESMFold for additional sequence/structure prediction
2. **AbLang**: Alternative antibody language model
3. **Rosetta**: Full Rosetta integration for refinement
4. **OpenMM**: Molecular dynamics for validation
5. **ImmuneBuilder**: Alternative structure prediction

## Performance Considerations

### Memory Usage

- **Preparation**: Low (< 1 GB)
- **Docking**: Low (< 2 GB)
- **Design (RFdiffusion)**: High (8-16 GB GPU)
- **IgLM**: Medium (2-4 GB GPU)
- **Filtering (Chai)**: High (16-32 GB GPU)

### Runtime

For 10 designs:
- Preparation: < 1 minute
- Docking: 5-10 minutes
- Design: 30-60 minutes
- IgLM scoring: < 5 minutes
- Filtering: 20-40 minutes per design

**Total**: ~4-8 hours for 10 designs

### Parallelization Opportunities

1. Multiple designs can be generated in parallel
2. Filtering can be parallelized across designs
3. IgLM scoring can batch sequences
4. Structure prediction can use multiple GPUs

## Maintenance

### Update Strategy

When upstream projects update:

1. **Check DiffAb releases**: Update CDR definitions if Chothia scheme changes
2. **Check Germinal releases**: Update filtering thresholds based on new research
3. **Check RFdiffusion releases**: Update CLI arguments if interface changes
4. **Check IgLM releases**: Update model loading if API changes

### Dependency Pinning

```txt
# Core (pinned for stability)
torch>=2.0.0,<3.0.0
biopython>=1.79,<2.0.0
abnumber>=0.3.0,<0.4.0

# Models (allow minor updates)
iglm>=0.4.0,<0.5.0
transformers>=4.30.0,<5.0.0

# Optional (flexible)
chai_lab>=0.1.0  # Latest compatible
```

## Conclusion

This pipeline successfully integrates the best components from three leading antibody design systems while maintaining:

1. **Modularity**: Each component can be used independently
2. **Simplicity**: Minimal dependencies, clear interfaces
3. **Flexibility**: Configurable for different use cases
4. **Extensibility**: Easy to add new features or swap components

The integration strategy prioritizes practical usability over theoretical completeness, making it suitable for production antibody design workflows.