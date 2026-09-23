#!/usr/bin/env bash
set -euo pipefail

WORKFLOW_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# These defaults define the confirmatory experiment. Values already exported
# by the caller remain authoritative, especially OUT_ROOT and cluster settings.
export LEAN_N="${LEAN_N:-12000}"
export LEAN_REPS="${LEAN_REPS:-100}"
export FACTOR_SPECS="${FACTOR_SPECS:-structure:5000}"
export PRIMARY_SE_METHOD="${PRIMARY_SE_METHOD:-structure}"
export PRIMARY_SE_RANK="${PRIMARY_SE_RANK:-5000}"
export LEAN_MAX_SIMULTANEOUS="${LEAN_MAX_SIMULTANEOUS:-30}"
export OUT_ROOT="${OUT_ROOT:-${WORKFLOW_DIR}/output_primary_r5000_100}"

exec bash "${WORKFLOW_DIR}/submit_predlmm_lean.sh" "$@"
