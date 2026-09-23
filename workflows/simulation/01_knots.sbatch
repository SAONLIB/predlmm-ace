#!/usr/bin/env bash
#SBATCH --job-name=sim_knots
#SBATCH --output=logs/01_knots_%A_%a.out
#SBATCH --error=logs/01_knots_%A_%a.err
set -euo pipefail; source "${WORKFLOW_DIR}/config.sh"; export PYTHONUNBUFFERED=1
row=$(awk -F'\t' -v id="$SLURM_ARRAY_TASK_ID" 'NR>1&&$1==id{print;exit}' "$OUT_ROOT/manifests/factors.tsv"); IFS=$'\t' read -r fid rid sim n h c rep frac rank method <<<"$row"
rdir=$(printf '%s/replicates/task_%06d' "$OUT_ROOT" "$rid"); fdir=$(printf '%s/factors/factor_%06d' "$OUT_ROOT" "$fid"); mkdir -p "$fdir"
"$PYTHON_EXE" "$WORKFLOW_DIR/select_knots.py" --metadata "$rdir/metadata.tsv" --rank "$rank" --method "$method" --seed "$((BASE_SEED+fid))" --out "$fdir/knot_ids.tsv"
