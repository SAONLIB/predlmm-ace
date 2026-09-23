#!/usr/bin/env bash
#SBATCH --job-name=sim_external
#SBATCH --output=logs/04_external_%A_%a.out
#SBATCH --error=logs/04_external_%A_%a.err
set -euo pipefail; source "${WORKFLOW_DIR}/config.sh"; export PYTHONUNBUFFERED=1
row=$(awk -F'\t' -v id="$SLURM_ARRAY_TASK_ID" 'NR>1&&$1==id{print;exit}' "$OUT_ROOT/manifests/replicates.tsv"); IFS=$'\t' read -r id sim n h c rep <<<"$row"
rdir=$(printf '%s/replicates/task_%06d' "$OUT_ROOT" "$id"); edir=$(printf '%s/external/task_%06d' "$OUT_ROOT" "$id")
"$PYTHON_EXE" "$WORKFLOW_DIR/external_methods.py" --data "$rdir" --grm-prefix "$rdir/grm/full" --simulation "$sim" --gcta "$GCTA_EXE" --threads "$THREADS" --cutoffs "$SPARSE_CUTOFFS" --out "$edir"
