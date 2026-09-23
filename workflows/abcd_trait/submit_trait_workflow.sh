#!/usr/bin/env bash
set -euo pipefail
WORKFLOW_DIR="$(cd "$(dirname "$0")" && pwd)"; export WORKFLOW_DIR
source "${WORKFLOW_DIR}/config_trait.sh"

for name in MASTER_METADATA BASE_OUT TRAIT_OUT GRM_PREFIX BFILE_PREFIX FACTOR_ROOT; do
  [[ -n "${!name:-}" ]] || { echo "Configuration value ${name} is required" >&2; exit 1; }
done
"${PYTHON_EXE}" -c 'import numpy, pandas, scipy, predlmm_ace' || {
  echo "Install this repository into PYTHON_EXE before submission" >&2
  exit 1
}
mkdir -p "${TRAIT_OUT}/data" "${TRAIT_OUT}/results" "${TRAIT_OUT}/logs"

for f in "${MASTER_METADATA}" "${GRM_PREFIX}.grm.id" "${BFILE_PREFIX}.bed" \
         "${FACTOR_ROOT}/rank4000_site_seed${SEED}/factor_metadata.json" \
         "${FACTOR_ROOT}/rank5000_site_seed${SEED}/factor_metadata.json"; do
  [[ -s "$f" ]] || { echo "Missing required input: $f" >&2; exit 1; }
done
"${PYTHON_EXE}" "${WORKFLOW_DIR}/prepare_trait.py" --metadata "${MASTER_METADATA}" \
  --grm-id "${GRM_PREFIX}.grm.id" --trait "${PHENOTYPE_COLUMN}" --family "${FAMILY_COLUMN}" \
  --continuous "${CONTINUOUS}" --categorical "${CATEGORICAL}" --out-dir "${TRAIT_OUT}"

if [[ "${1:-}" == "--dry-run" ]]; then
  echo "Validated and prepared ${TRAIT}; would submit 4 PredLMM tasks, 1 GCTA task, 2 fastGWA tasks, and summary."
  exit 0
fi

pred=$(sbatch --parsable --export=ALL,WORKFLOW_DIR="${WORKFLOW_DIR}" --array=0-3 \
  --cpus-per-task="${PRED_CPUS}" --mem="${PRED_MEM}" --time="${PRED_TIME}" \
  --output="${TRAIT_OUT}/logs/predlmm_%A_%a.out" --error="${TRAIT_OUT}/logs/predlmm_%A_%a.err" \
  "${WORKFLOW_DIR}/01_predlmm.sbatch")
gcta=$(sbatch --parsable --export=ALL,WORKFLOW_DIR="${WORKFLOW_DIR}" \
  --cpus-per-task="${GCTA_CPUS}" --mem="${GCTA_MEM}" --time="${GCTA_TIME}" \
  --output="${TRAIT_OUT}/logs/gcta_%j.out" --error="${TRAIT_OUT}/logs/gcta_%j.err" \
  "${WORKFLOW_DIR}/02_gcta.sbatch")
nfast=$(wc -w <<<"${FASTGWA_LABELS}"); last=$((nfast-1))
fast=$(sbatch --parsable --export=ALL,WORKFLOW_DIR="${WORKFLOW_DIR}" --array="0-${last}" \
  --cpus-per-task="${FASTGWA_CPUS}" --mem="${FASTGWA_MEM}" --time="${FASTGWA_TIME}" \
  --output="${TRAIT_OUT}/logs/fastgwa_%A_%a.out" --error="${TRAIT_OUT}/logs/fastgwa_%A_%a.err" \
  "${WORKFLOW_DIR}/03_fastgwa.sbatch")
summary=$(sbatch --parsable --export=ALL,WORKFLOW_DIR="${WORKFLOW_DIR}",PRED_JOB="${pred}",GCTA_JOB="${gcta}",FAST_JOB="${fast}" \
  --dependency="afterany:${pred}:${gcta}:${fast}" --output="${TRAIT_OUT}/logs/summary_%j.out" \
  --error="${TRAIT_OUT}/logs/summary_%j.err" "${WORKFLOW_DIR}/04_summary.sbatch")
printf 'stage\tjob_id\nPredLMM\t%s\nGCTA\t%s\nfastGWA\t%s\nsummary\t%s\n' "$pred" "$gcta" "$fast" "$summary" | tee "${TRAIT_OUT}/workflow_jobs.tsv"
