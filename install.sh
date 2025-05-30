#!/bin/bash

set -e

echo "==========================================="
echo "Phytovenomics Installation Script"
echo "==========================================="

# Create necessary directories
echo "Creating directories..."
mkdir -p config models data/temp logs results

# Check if Python 3 is installed
if command -v python3 &>/dev/null; then
    PYTHON="python3"
elif command -v python &>/dev/null; then
    # Check if this Python is version 3+
    PY_VERSION=$(python --version 2>&1)
    if [[ $PY_VERSION == *"Python 3"* ]]; then
        PYTHON="python"
    else
        echo "Error: Python 3 is required but not found"
        exit 1
    fi
else
    echo "Error: Python is not installed"
    exit 1
fi

echo "Using Python: $($PYTHON --version)"

# Check if pip is installed
if command -v pip3 &>/dev/null; then
    PIP="pip3"
elif command -v pip &>/dev/null; then
    PIP="pip"
else
    echo "Error: pip is not installed"
    exit 1
fi

echo "Using Pip: $($PIP --version)"

# Install the package
echo "\nInstalling Phytovenomics..."
$PIP install -e .

# Initialize configuration
echo "\nInitializing configuration..."
$PYTHON -m phytovenomics.cli.setup_cli init_config

# Run setup
echo "\nWould you like to run the full setup now? This includes downloading models and installing dependencies. [y/N]"
read -r run_setup

if [[ $run_setup =~ ^[Yy]$ ]]; then
    echo "Running setup..."
    $PYTHON -m phytovenomics.cli.setup_cli setup
else
    echo "\nSkipping full setup."
    echo "You can run the setup later with: python -m phytovenomics.cli.setup_cli setup"
fi

echo "\n==========================================="
echo "Installation complete!"
echo "==========================================="
echo "To get started, run: phytovenomics-run demo"
echo "For more information, run: phytovenomics-setup --help"
echo "=========================================="
