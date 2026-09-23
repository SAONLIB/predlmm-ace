#!/usr/bin/env bash
#SBATCH --job-name=sim_factor
#SBATCH --output=logs/02_factor_%A_%a.out
#SBATCH --error=logs/02_factor_%A_%a.err
set -euo pipefail; source "${WORKFLOW_DIR}/config.sh"; export PYTHONUNBUFFERED=1
row=$(awk -F'\t' -v id="$SLURM_ARRAY_TASK_ID" 'NR>1&&$1==id{print;exit}' "$OUT_ROOT/manifests/factors.tsv"); IFS=$'\t' read -r fid rid sim n h c rep frac rank method <<<"$row"
rdir=$(printf '%s/replicates/task_%06d' "$OUT_ROOT" "$rid"); fdir=$(printf '%s/factors/factor_%06d' "$OUT_ROOT" "$fid")
"$PYTHON_EXE" -m predlmm_ace.grm_to_nystrom --grm-prefix "$rdir/grm/full" --knots "$fdir/knot_ids.tsv" --out "$fdir"
