#!/usr/bin/env bash
set -euo pipefail
WORKFLOW_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export WORKFLOW_DIR
requested_out_root="${OUT_ROOT:-}"
source "${WORKFLOW_DIR}/config.sh"
[[ -n "${requested_out_root}" ]] || { echo "Set OUT_ROOT first" >&2; exit 1; }
export OUT_ROOT="${requested_out_root}"
"${PYTHON_EXE}" "${WORKFLOW_DIR}/summarize_lean.py" --root "${OUT_ROOT}"
