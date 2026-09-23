# Simulation workflows

Install the repository first (`python -m pip install -e .`) and make GCTA
available. Copy `config.example.sh` to `config.local.sh` and set cluster paths.

## Confirmatory primary experiment

```bash
bash submit_primary_100.sh --dry-run
bash submit_primary_100.sh
```

This submits six scenarios at `N=12,000`, 100 replicates per scenario, and only
the structure-aware rank-5,000 PredLMM-ACE configuration. All 600 fits request
profile-Hessian standard errors. The dry run enforces the expected 600/600 task
counts.

## Rank and knot sensitivity

```bash
export OUT_ROOT=/path/to/a/new/sensitivity/output
export LEAN_REPS=20
export FACTOR_SPECS='structure:3000,structure:4000,structure:5000,random:5000'
bash submit_predlmm_lean.sh --dry-run
bash submit_predlmm_lean.sh
```

Only structure-aware rank 5,000 receives Hessian SEs by default. Change
`PRIMARY_SE_METHOD` and `PRIMARY_SE_RANK` only when another primary
configuration is intended.

`submit_workflow.sh` is the expanded workflow that also schedules external
GCTA/fastGWA comparators. These jobs can be much slower and are not required
for the 100-replicate primary confirmation.
