#!/bin/bash

# setup_igfold.sh: Script to set up IgFold dependencies and run tests

echo "Setting up IgFold integration environment..."

# Create directories if they don't exist
mkdir -p src/structure_prediction/results

# Install dependencies
echo "Installing dependencies from requirements_igfold.txt..."
pip install -r requirements_igfold.txt

# Run basic availability test
echo "Testing model availability..."
python -c "from src.structure_prediction.hybrid_prediction import IgFoldInterface, ESMFoldInterface; \
print('IgFold available:', IgFoldInterface().is_available()); \
print('ESMFold available:', ESMFoldInterface().is_available())"

# Run the test script
echo "Running hybrid prediction tests..."
python src/structure_prediction/test_hybrid_prediction.py --quick

echo "Setup complete."
