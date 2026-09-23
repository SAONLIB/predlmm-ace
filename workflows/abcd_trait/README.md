# Reusable ABCD trait workflow

This workflow reuses an authorized cohort-wide GRM, PLINK bed set, sparse GRMs,
and precomputed rank-4,000/rank-5,000 Nyström factors. It prepares a complete-
case trait sample, selects one deterministic participant per family, fits
PredLMM-ACE/AE with observed-Hessian SEs, runs sample-matched GCTA-GREML and
fastGWA sensitivity analyses, and creates CSV/LaTeX comparison tables.

No ABCD data, participant IDs, results, or institution-specific paths are
included in this repository.

## Configure

```bash
cp config.local.example.sh config.local.sh
```

Set at least:

```bash
export TRAIT="your_trait"
export MASTER_METADATA="/authorized/path/metadata.tsv"
export BASE_OUT="/authorized/project/output"
export GRM_PREFIX="/authorized/path/full_grm_prefix"
export BFILE_PREFIX="/authorized/path/plink_prefix"
export FACTOR_ROOT="/authorized/path/nystrom"
export GCTA="/path/to/gcta64"
export PYTHON_EXE="/path/to/python"
```

The metadata must contain unique `FID`, `IID`, the phenotype, `family_id`,
`age_years`, `sex`, `site_id`, and `PC1`--`PC10`. Update `CONTINUOUS` and
`CATEGORICAL` when the analysis plan differs.

## Validate and submit

```bash
bash submit_trait_workflow.sh --dry-run
bash submit_trait_workflow.sh
```

The dry run prepares and validates inputs but submits no jobs. Job IDs are
written to `${TRAIT_OUT}/workflow_jobs.tsv`; final tables are placed in
`${TRAIT_OUT}/results/`.

The SEs condition on the chosen rank and knots. The full-sample AE model is a
misspecification diagnostic; the one-participant-per-family AE analysis is the
sample-matched comparison with GCTA.
