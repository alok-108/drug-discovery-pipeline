#!/usr/bin/env bash
# ==============================================================================
# Script: setup_env.sh
# Purpose: Initialize Conda/Python virtual environment and install dependencies
# ==============================================================================

set -euo pipefail

echo "======================================================================"
echo " Setting up AI Drug Discovery Pipeline Environment..."
echo "======================================================================"

if command -v conda &> /dev/null; then
    echo "Conda detected. Creating environment from environment.yml..."
    conda env create -f environment.yml || conda env update -f environment.yml
    echo ""
    echo "To activate environment, run:"
    echo "  conda activate drug-discovery"
else
    echo "Conda not found. Creating standard Python virtual environment..."
    python3 -m venv .venv
    source .venv/bin/activate
    pip install --upgrade pip
    pip install -r requirements.txt
    pip install -e .
    echo ""
    echo "To activate virtual environment, run:"
    echo "  source .venv/bin/activate"
fi

echo "======================================================================"
echo " Environment setup complete! Run 'pytest tests/ -v' to verify."
echo "======================================================================"
