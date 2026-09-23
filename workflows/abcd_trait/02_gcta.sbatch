#!/usr/bin/env bash
set -euo pipefail
source "${WORKFLOW_DIR}/config_trait.sh"
KEEP="${TRAIT_OUT}/data/${TRAIT}_onefamily.keep"
OUT="${TRAIT_OUT}/results/${TRAIT}_GCTA_onefamily"
"${GCTA}" --grm "${GRM_PREFIX}" --keep "${KEEP}" \
  --pheno "${TRAIT_OUT}/data/${TRAIT}.pheno" \
  --qcovar "${TRAIT_OUT}/data/${TRAIT}.qcovar" \
  --covar "${TRAIT_OUT}/data/${TRAIT}.covar" \
  --reml --thread-num "${GCTA_CPUS}" --out "${OUT}"
test -s "${OUT}.hsq"

