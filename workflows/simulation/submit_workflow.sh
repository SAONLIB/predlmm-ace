#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/config.sh"
dry=0; [[ "${1:-}" == '--dry-run' ]] && dry=1
mkdir -p "$OUT_ROOT"/{logs,manifests,replicates,factors,results,external}
"$PYTHON_EXE" "$WORKFLOW_DIR/make_manifests.py" --reps "$REPS" --fractions "$KNOT_FRACTIONS" --methods "$KNOT_METHODS" --out "$OUT_ROOT/manifests"
nrep=$(($(wc -l < "$OUT_ROOT/manifests/replicates.tsv")-1)); nfac=$(($(wc -l < "$OUT_ROOT/manifests/factors.tsv")-1))
echo "replicate tasks=$nrep factor/fit tasks=$nfac output=$OUT_ROOT"
((dry)) && { echo 'Dry run only'; exit 0; }
opts=(--parsable --cpus-per-task="$SLURM_CPUS" --mem="$SLURM_MEM" --time="$SLURM_TIME")
[[ -n "$SLURM_ACCOUNT" ]] && opts+=(--account="$SLURM_ACCOUNT")
[[ -n "$SLURM_PARTITION" ]] && opts+=(--partition="$SLURM_PARTITION")
submit_chunks(){ local count=$1 script=$2 dep=${3:-} start=0 ids=() end jid arr; while ((start<count)); do end=$((start+899)); ((end>=count))&&end=$((count-1)); arr="${start}-${end}%${MAX_SIMULTANEOUS}"; args=("${opts[@]}" --array="$arr"); [[ -n "$dep" ]]&&args+=(--dependency="afterok:${dep}"); jid=$(cd "$OUT_ROOT"&&sbatch "${args[@]}" "$script"); ids+=("$jid"); start=$((end+1)); done; (IFS=:; echo "${ids[*]}"); }
sim=$(submit_chunks "$nrep" "$WORKFLOW_DIR/00_simulate.sbatch")
grm=$(submit_chunks "$nrep" "$WORKFLOW_DIR/01_make_grm.sbatch" "$sim")
knots=$(submit_chunks "$nfac" "$WORKFLOW_DIR/01_knots.sbatch" "$sim")
factor=$(submit_chunks "$nfac" "$WORKFLOW_DIR/02_factor.sbatch" "${grm}:${knots}")
fit=$(submit_chunks "$nfac" "$WORKFLOW_DIR/03_fit_se.sbatch" "$factor")
external=$(submit_chunks "$nrep" "$WORKFLOW_DIR/04_external.sbatch" "$grm")
summary=$(cd "$OUT_ROOT"&&sbatch --parsable --dependency="afterok:${fit}:${external}" "$WORKFLOW_DIR/05_summary.sbatch")
printf 'stage\tjob_ids\nsimulate\t%s\ngrm\t%s\nknots\t%s\nfactor\t%s\nfit_se\t%s\nexternal\t%s\nsummary\t%s\n' "$sim" "$grm" "$knots" "$factor" "$fit" "$external" "$summary" | tee "$OUT_ROOT/workflow_jobs.tsv"

