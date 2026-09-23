#!/usr/bin/env bash
#SBATCH --job-name=sim_generate
#SBATCH --output=logs/00_sim_%A_%a.out
#SBATCH --error=logs/00_sim_%A_%a.err
set -euo pipefail; source "${WORKFLOW_DIR}/config.sh"; export PYTHONUNBUFFERED=1
row=$(awk -F'\t' -v id="$SLURM_ARRAY_TASK_ID" 'NR>1&&$1==id{print;exit}' "$OUT_ROOT/manifests/replicates.tsv"); IFS=$'\t' read -r id sim n h c rep <<<"$row"
[[ -n "$row" ]] || { echo "No replicate manifest row for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
h_code=$(awk -v x="$h" 'BEGIN {printf "%d", 100*x+0.5}')
c_code=$(awk -v x="$c" 'BEGIN {printf "%d", 100*x+0.5}')
seed=$((BASE_SEED + sim*100000000 + n*1000 + h_code*10000 + c_code*100 + rep))
dir=$(printf '%s/replicates/task_%06d' "$OUT_ROOT" "$id"); mkdir -p "$dir"
echo "task=$id simulation=$sim N=$n h2=$h c2=$c replicate=$rep seed=$seed"
"$PYTHON_EXE" "$WORKFLOW_DIR/simulate_replicate.py" --simulation "$sim" --N "$n" --M "$M" --h2 "$h" --c2 "$c" --replicate "$rep" --seed "$seed" --n-pcs "$N_PCS" --out "$dir"
