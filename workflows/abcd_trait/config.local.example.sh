#!/usr/bin/env bash
# Copy to config.local.sh. Do not commit authorized-data paths.

export TRAIT="your_trait"
export MASTER_METADATA="/authorized/path/metadata.tsv"
export BASE_OUT="/authorized/project/output"
export TRAIT_OUT="${BASE_OUT}/trait_analyses/${TRAIT}"

export GRM_PREFIX="/authorized/path/full_grm_prefix"
export BFILE_PREFIX="/authorized/path/plink_prefix"
export FACTOR_ROOT="/authorized/path/nystrom"

export FASTGWA_PREFIX_0050="/authorized/path/sparse_0050"
export FASTGWA_PREFIX_0100="/authorized/path/sparse_0100"

export GCTA="/path/to/gcta64"
export PYTHON_EXE="/path/to/python"
