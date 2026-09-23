#!/usr/bin/env bash
set -euo pipefail
source "${WORKFLOW_DIR}/config_trait.sh"
if [[ -n "${PRED_JOB:-}" && -n "${GCTA_JOB:-}" && -n "${FAST_JOB:-}" ]]; then
  sacct -j "${PRED_JOB},${GCTA_JOB},${FAST_JOB}" --units=G \
    --format=JobID%24,JobName%22,State%15,ExitCode,Elapsed,TotalCPU,AllocCPUS,ReqMem,MaxRSS \
    > "${TRAIT_OUT}/results/${TRAIT}_computational_efficiency.txt" || true
fi
"${PYTHON_EXE}" "${WORKFLOW_DIR}/summarize_trait.py" --trait "${TRAIT}" --out-dir "${TRAIT_OUT}"
