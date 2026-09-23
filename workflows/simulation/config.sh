#!/usr/bin/env bash
set -euo pipefail
export WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export OUT_ROOT="${OUT_ROOT:-${WORKFLOW_DIR}/output}"
export PYTHON_EXE="${PYTHON_EXE:-python}"
export GCTA_EXE="${GCTA_EXE:-gcta64}"
export REPS="${REPS:-100}"
export M="${M:-16000}"
export N_PCS="${N_PCS:-5}"
export KNOT_FRACTIONS="${KNOT_FRACTIONS:-0.10,0.25}"
export KNOT_METHODS="${KNOT_METHODS:-structure,random}"
export SPARSE_CUTOFFS="${SPARSE_CUTOFFS:-0.025,0.05,0.10}"
export BASE_SEED="${BASE_SEED:-20260921}"
export THREADS="${THREADS:-8}"
export MAX_SIMULTANEOUS="${MAX_SIMULTANEOUS:-4}"
export SLURM_TIME="${SLURM_TIME:-24:00:00}"
export SLURM_MEM="${SLURM_MEM:-64G}"
export SLURM_CPUS="${SLURM_CPUS:-8}"
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"
export SLURM_PARTITION="${SLURM_PARTITION:-}"

# Optional site-specific settings. This file is ignored by Git.
if [[ -s "${WORKFLOW_DIR}/config.local.sh" ]]; then
    # shellcheck source=/dev/null
    source "${WORKFLOW_DIR}/config.local.sh"
fi
