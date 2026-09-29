#!/usr/bin/env bash
# ==============================================================================
# Script: run_pipeline.sh
# Purpose: Execute end-to-end AI drug discovery pipeline
# ==============================================================================

set -euo pipefail

echo "======================================================================"
echo " Starting AI-Driven Drug Discovery Pipeline..."
echo "======================================================================"

# Determine project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

# Set PYTHONPATH
export PYTHONPATH="${PROJECT_ROOT}/src:${PYTHONPATH:-}"

CONFIG_FILE="${1:-configs/pipeline_config.yaml}"

echo "Using configuration: $CONFIG_FILE"
echo "Running full pipeline orchestrator..."

python3 -c "from utils.io_utils import run_full_pipeline_cli; run_full_pipeline_cli('$CONFIG_FILE')"

echo "======================================================================"
echo " Pipeline execution completed successfully."
echo " Check outputs/ directory for prioritized targets, pockets, and designs."
echo "======================================================================"
