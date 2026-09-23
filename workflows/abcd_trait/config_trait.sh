#!/usr/bin/env bash
# Put site-specific values in config.local.sh. Restricted paths and data must
# never be committed to Git.

CONFIG_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if [[ -s "${CONFIG_DIR}/config.local.sh" ]]; then
  # shellcheck source=/dev/null
  source "${CONFIG_DIR}/config.local.sh"
fi

export TRAIT="${TRAIT:-example_trait}"
export SEED="${SEED:-20260918}"

# Master metadata must contain FID, IID, the phenotype, family_id, age_years,
# sex, site_id, and PC1-PC10. Rows need not be in GRM order.
export MASTER_METADATA="${MASTER_METADATA:-}"

export BASE_OUT="${BASE_OUT:-}"
export TRAIT_OUT="${TRAIT_OUT:-${BASE_OUT:+${BASE_OUT}/trait_analyses/${TRAIT}}}"

# Reused genotype resources. These must contain every candidate participant.
export GRM_PREFIX="${GRM_PREFIX:-${BASE_OUT:+${BASE_OUT}/grm/full}}"
export BFILE_PREFIX="${BFILE_PREFIX:-${BASE_OUT:+${BASE_OUT}/plink/full}}"
export FACTOR_ROOT="${FACTOR_ROOT:-${BASE_OUT:+${BASE_OUT}/nystrom}}"

# Existing sparse GRMs for fastGWA sensitivity analyses.
export FASTGWA_LABELS="${FASTGWA_LABELS:-0050 0100}"
export FASTGWA_PREFIX_0050="${FASTGWA_PREFIX_0050:-${BASE_OUT:+${BASE_OUT}/grm/sparse_0050}}"
export FASTGWA_PREFIX_0100="${FASTGWA_PREFIX_0100:-${BASE_OUT:+${BASE_OUT}/grm/sparse_0100}}"

export GCTA="${GCTA:-gcta64}"
export PYTHON_EXE="${PYTHON_EXE:-python}"

export PHENOTYPE_COLUMN="${PHENOTYPE_COLUMN:-${TRAIT}}"
export FAMILY_COLUMN="${FAMILY_COLUMN:-family_id}"
export CONTINUOUS="${CONTINUOUS:-age_years,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10}"
export CATEGORICAL="${CATEGORICAL:-sex,site_id}"

export PRED_CPUS="${PRED_CPUS:-4}"
export PRED_MEM="${PRED_MEM:-64G}"
export PRED_TIME="${PRED_TIME:-06:00:00}"
export GCTA_CPUS="${GCTA_CPUS:-16}"
export GCTA_MEM="${GCTA_MEM:-48G}"
export GCTA_TIME="${GCTA_TIME:-06:00:00}"
export FASTGWA_CPUS="${FASTGWA_CPUS:-8}"
export FASTGWA_MEM="${FASTGWA_MEM:-24G}"
export FASTGWA_TIME="${FASTGWA_TIME:-12:00:00}"
