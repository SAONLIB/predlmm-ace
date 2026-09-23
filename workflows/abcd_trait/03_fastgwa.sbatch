#!/usr/bin/env bash
set -euo pipefail
source "${WORKFLOW_DIR}/config_trait.sh"
read -r -a LABELS <<<"${FASTGWA_LABELS}"
LABEL="${LABELS[$SLURM_ARRAY_TASK_ID]}"
case "${LABEL}" in
  0050) SPARSE="${FASTGWA_PREFIX_0050}" ;;
  0100) SPARSE="${FASTGWA_PREFIX_0100}" ;;
  *) echo "Unsupported sparse-GRM label: ${LABEL}" >&2; exit 2 ;;
esac
OUT="${TRAIT_OUT}/results/${TRAIT}_fastGWA_${LABEL}"
"${GCTA}" --bfile "${BFILE_PREFIX}" --grm-sparse "${SPARSE}" \
  --keep "${TRAIT_OUT}/data/${TRAIT}_full.keep" \
  --pheno "${TRAIT_OUT}/data/${TRAIT}.pheno" \
  --qcovar "${TRAIT_OUT}/data/${TRAIT}.qcovar" \
  --covar "${TRAIT_OUT}/data/${TRAIT}.covar" \
  --fastGWA-mlm --thread-num "${FASTGWA_CPUS}" --out "${OUT}"
test -s "${OUT}.log"

