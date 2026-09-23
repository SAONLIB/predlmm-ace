#!/usr/bin/env bash
# Four tasks: ACE r4000/r5000, full-sample AE r5000, one-family AE r5000.
set -euo pipefail
source "${WORKFLOW_DIR}/config_trait.sh"
mkdir -p "${TRAIT_OUT}/results" "${TRAIT_OUT}/logs"

case "${SLURM_ARRAY_TASK_ID}" in
  0) MODEL=ACE; RANK=4000; SAMPLE=full ;;
  1) MODEL=ACE; RANK=5000; SAMPLE=full ;;
  2) MODEL=AE;  RANK=5000; SAMPLE=full ;;
  3) MODEL=AE;  RANK=5000; SAMPLE=onefamily ;;
  *) echo "Invalid task" >&2; exit 2 ;;
esac

FACTOR_DIR="${FACTOR_ROOT}/rank${RANK}_site_seed${SEED}"
METADATA="${TRAIT_OUT}/data/${TRAIT}_metadata.tsv"
KEEP="${TRAIT_OUT}/data/${TRAIT}_${SAMPLE}.keep"
BASE="${TRAIT_OUT}/results/${TRAIT}_PredLMM_${MODEL}_${SAMPLE}_rank${RANK}"

for f in "${FACTOR_DIR}/factor_metadata.json" "${METADATA}" "${KEEP}"; do
  [[ -s "$f" ]] || { echo "Missing $f" >&2; exit 1; }
done

"${PYTHON_EXE}" -m predlmm_ace.fit \
  --factor-dir "${FACTOR_DIR}" --metadata "${METADATA}" \
  --phenotype "${PHENOTYPE_COLUMN}" --family-column "${FAMILY_COLUMN}" \
  --continuous "${CONTINUOUS}" --categorical "${CATEGORICAL}" \
  --keep "${KEEP}" --model "${MODEL}" --out "${BASE}.json"

"${PYTHON_EXE}" -m predlmm_ace.profile_se \
  --fit-json "${BASE}.json" --factor-dir "${FACTOR_DIR}" \
  --metadata "${METADATA}" --phenotype "${PHENOTYPE_COLUMN}" \
  --family-column "${FAMILY_COLUMN}" --continuous "${CONTINUOUS}" \
  --categorical "${CATEGORICAL}" --keep "${KEEP}" \
  --hessian-step 1e-4 --out "${BASE}_se.json"
