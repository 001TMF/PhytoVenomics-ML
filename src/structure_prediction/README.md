# Hybrid Structure Prediction Module

This module implements a hybrid approach to protein structure prediction, specifically optimized for antibody structures. By combining the strengths of general protein structure prediction models (ESMFold) with specialized antibody structure prediction models (IgFold), we achieve higher accuracy predictions especially for CDR regions.

## 🧬 Overview

The hybrid prediction pipeline:
1. Offers specialized predictions for antibody structures
2. Falls back to general protein structure prediction when specialized models are unavailable
3. Intelligently selects the best model based on sequence type and available resources
4. Caches predictions to improve performance for repeated queries
5. Provides quality assessment and visualization tools

## 📋 Requirements

### Core Dependencies
- Python 3.8+
- PyTorch 1.12+
- NumPy
- BioPython
- Matplotlib

### Optional Dependencies
- ESMFold for general protein structure prediction
- IgFold for specialized antibody structure prediction
- py3Dmol for interactive 3D visualization

The module is designed with graceful fallbacks when dependencies are missing:
- Without ESMFold: General structure prediction will be unavailable
- Without IgFold: Specialized antibody prediction will be unavailable 
- Without py3Dmol: Visualization will default to static Matplotlib plots

## 🚀 Installation

### Basic Installation

```bash
# Clone the repository
git clone https://github.com/your-organization/phytovenomics.git
cd phytovenomics

# Install basic dependencies
pip install -r requirements.txt
```

### Full Installation (with ESMFold and IgFold)

```bash
# Install ESMFold
pip install "fair-esm[esmfold]"

# Install IgFold
git clone https://github.com/Graylab/IgFold.git
cd IgFold
pip install -e .
cd ..

# Install py3Dmol for interactive visualization
pip install py3Dmol
```

Alternatively, you can use our setup script:

```bash
bash setup_igfold.sh
```

## 🏄‍♂️ Usage Examples

### Basic Structure Prediction

```python
from structure_prediction.hybrid_prediction import HybridPredictionPipeline

# Initialize the pipeline
pipeline = HybridPredictionPipeline(use_cache=True)

# Predict antibody structure
heavy_chain_sequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKGWFDYWGQGTLVTVSS"

result = pipeline.predict_antibody_structure(
    sequence=heavy_chain_sequence,
    chain_type="heavy",  # Options: "heavy", "light", or "paired"
    use_ensemble=True    # Use both models when available
)

# Analyze structure quality
quality_metrics = pipeline.analyze_structure_quality(result)
print(f"Confidence: {result.confidence:.4f}")

# Visualize structure
viewer = pipeline.visualize_structure(
    result, 
    highlight_cdrs=True,
    visualization_type="py3Dmol"  # or "matplotlib" for static image
)
display(viewer)
```

### Paired Antibody Chain Prediction

```python
from structure_prediction.hybrid_prediction import HybridPredictionPipeline

# Initialize the pipeline
pipeline = HybridPredictionPipeline()

# Predict paired antibody structure (heavy and light chain together)
heavy_chain = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS..."
light_chain = "DIVMTQSPLSLPVTPGEPASISCRSSQSLL..."
paired_sequence = heavy_chain + light_chain  # Simple concatenation for demonstration

result = pipeline.predict_antibody_structure(
    sequence=paired_sequence,
    chain_type="paired",
    prefer_model="igfold"  # Prefer IgFold for paired prediction
)

# Save visualization to file
pipeline.visualize_structure(
    result,
    output_path="antibody_structure.png",
    highlight_cdrs=True
)
print(f"Structure saved to antibody_structure.png")
```

### Handling Missing Dependencies

```python
from structure_prediction.hybrid_prediction import HybridPredictionPipeline

# Initialize the pipeline - will work with whatever models are available
try:
    pipeline = HybridPredictionPipeline()
    
    # Check which models are available
    if pipeline.esm_interface.is_available():
        print("ESMFold is available for general protein structure prediction")
    if pipeline.igfold_interface.is_available():
        print("IgFold is available for specialized antibody structure prediction")
    
    # Proceed with prediction if any model is available
    sequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS..."
    result = pipeline.predict_antibody_structure(sequence)
    
 except RuntimeError as e:
    print(f"Structure prediction unavailable: {e}")
    print("Please install ESMFold or IgFold to enable structure prediction")
```

## 📊 Performance Considerations

- **GPU Acceleration**: Both ESMFold and IgFold benefit significantly from GPU acceleration. For optimal performance, we recommend using a CUDA-compatible GPU.

- **Memory Requirements**: Large protein complexes may require substantial memory. ESMFold typically requires 16+ GB of GPU memory for full-sized proteins.

- **Caching**: The pipeline automatically caches predictions to improve performance for repeated queries. To disable caching, initialize the pipeline with `use_cache=False`.

- **Model Selection**: For antibody structures, IgFold typically produces higher quality predictions, especially for CDR regions. For other proteins, ESMFold is recommended.

## 🔍 Troubleshooting

### Common Issues

1. **ImportError for ESM or IgFold**
   - Ensure you've installed all dependencies according to the installation instructions
   - Check CUDA compatibility for PyTorch if using GPU acceleration

2. **CUDA Out of Memory**
   - Try reducing batch size or using a smaller model variant
   - Process one chain at a time instead of paired prediction

3. **Visualization Issues**
   - If py3Dmol is not available, the pipeline will automatically fall back to Matplotlib
   - Check that your Jupyter notebook environment supports interactive visualization

## 📝 Citation

If you use this module in your research, please cite:

- ESMFold: Lin, Z. et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379, 1123–1130 (2023).
- IgFold: Ruffolo, J.A., et al. Fast, accurate antibody structure prediction from deep learning on massive set of natural antibodies. Nature Methods (2023).

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.
