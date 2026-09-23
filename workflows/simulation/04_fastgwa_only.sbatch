#!/usr/bin/env bash
#SBATCH --job-name=sim_fastgwa
#SBATCH --output=logs/04_fastgwa_%A_%a.out
#SBATCH --error=logs/04_fastgwa_%A_%a.err

set -euo pipefail
source "${WORKFLOW_DIR}/config.sh"
export PYTHONUNBUFFERED=1
FASTGWA_CUTOFF="${FASTGWA_CUTOFF:-0.05}"

row=$(awk -F'\t' -v id="${SLURM_ARRAY_TASK_ID}" \
    'NR > 1 && $1 == id {print; exit}' \
    "${OUT_ROOT}/manifests/replicates.tsv")
IFS=$'\t' read -r id sim n h c rep <<< "${row}"
rdir=$(printf '%s/replicates/task_%06d' "${OUT_ROOT}" "${id}")
edir=$(printf '%s/external/task_%06d' "${OUT_ROOT}" "${id}")

"${PYTHON_EXE}" "${WORKFLOW_DIR}/fastgwa_only.py" \
    --data "${rdir}" \
    --grm-prefix "${rdir}/grm/full" \
    --gcta "${GCTA_EXE}" \
    --threads "${THREADS}" \
    --cutoff "${FASTGWA_CUTOFF}" \
    --out "${edir}"
