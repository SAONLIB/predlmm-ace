#!/usr/bin/env bash
#SBATCH --job-name=sim_summary
#SBATCH --output=logs/05_summary_%j.out
#SBATCH --error=logs/05_summary_%j.err

set -euo pipefail
source "${WORKFLOW_DIR}/config.sh"
"${PYTHON_EXE}" "${WORKFLOW_DIR}/summarize_lean.py" --root "${OUT_ROOT}"
