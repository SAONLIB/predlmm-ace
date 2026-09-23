#!/usr/bin/env bash
#SBATCH --job-name=sim_grm
#SBATCH --output=logs/01_grm_%A_%a.out
#SBATCH --error=logs/01_grm_%A_%a.err
set -euo pipefail
source "${WORKFLOW_DIR}/config.sh"
export PYTHONUNBUFFERED=1
row=$(awk -F'\t' -v id="$SLURM_ARRAY_TASK_ID" 'NR>1&&$1==id{print;exit}' "$OUT_ROOT/manifests/replicates.tsv")
IFS=$'\t' read -r id sim n h c rep <<<"$row"
rdir=$(printf '%s/replicates/task_%06d' "$OUT_ROOT" "$id")
mkdir -p "$rdir/grm"
"$GCTA_EXE" --bfile "$rdir/geno" --make-grm --thread-num "$THREADS" --out "$rdir/grm/full"
echo 'COMPLETE GRM'
