# PredLMM-ACE

PredLMM-ACE estimates additive-genetic, shared-family or household, and
residual variance components using a low-rank Nyström representation of a
genetic relationship matrix (GRM):

\[
V = \sigma_g^2\{BB^T + \operatorname{diag}(d)\}
  + \sigma_c^2\Pi + \sigma_e^2 I.
\]

The exact backend profiles the overall variance and uses Woodbury identities;
it never forms an `N x N` covariance or inverse. Family covariance is stored as
a label vector and evaluated analytically by family block. The release also
contains the simulation and ABCD-trait Slurm workflows used for the manuscript.

## Status

- `ExactProfileML` is the validated backend used for the reported simulations
  and ABCD analysis.
- Observed profile-Hessian standard errors for `h2` and `c2` are provided for
  interior solutions.
- `MatrixFreeAIREML` is experimental and is not the estimator used for the
  reported results. It requires additional large-cohort validation before
  definitive inference.
- Ordinary Wald inference is not valid when a variance component is on the
  boundary. Retain the point estimate and use a boundary-aware procedure.

## Installation

```bash
cd predlmm-ace
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[test]'
pytest
```

GCTA is required only for workflows that construct a GRM or run GCTA/fastGWA
comparators. The core Python package does not require GCTA.

## Command-line tools

```bash
predlmm-select-knots --help
predlmm-grm-to-nystrom --help
predlmm-fit --help
predlmm-profile-se --help
```

The input metadata must contain unique `FID`/`IID` pairs and the requested
phenotype, family label, and covariates. Factor sample order is authoritative;
all metadata are explicitly aligned to `sample_ids.tsv` before fitting.

## Reproduce the manuscript simulations

The recommended final design separates the confirmatory and sensitivity runs:

```bash
cd workflows/simulation
cp config.example.sh config.local.sh
# Edit paths and Slurm settings in config.local.sh, then source it.
source config.local.sh

# Primary analysis: six scenarios, 100 replicates/scenario,
# structure-aware rank 5,000, with SEs for every fit.
bash submit_primary_100.sh --dry-run
bash submit_primary_100.sh

# Optional 20-replicate rank/knot sensitivity analysis.
export LEAN_REPS=20
export FACTOR_SPECS='structure:3000,structure:4000,structure:5000,random:5000'
export OUT_ROOT=/path/to/sensitivity_output
bash submit_predlmm_lean.sh --dry-run
bash submit_predlmm_lean.sh
```

The primary run creates 600 simulated data sets and 600 fits. The sensitivity
run creates 120 data sets and 480 fits. External GCTA/fastGWA jobs are excluded
from the lean workflow because they are substantially slower; separate scripts
are included for targeted benchmarking.

## Reusable ABCD trait workflow

See [`workflows/abcd_trait/README.md`](workflows/abcd_trait/README.md). Restricted
ABCD data, paths, IDs, and results are not included. The workflow expects an
authorized local metadata table and existing genotype/GRM resources.

## Repository layout

```text
src/predlmm_ace/       statistical implementation and command-line tools
tests/                 dense-algebra and Nyström equivalence tests
workflows/simulation/  Slurm simulation workflows and summaries
workflows/abcd_trait/  reusable authorized-data trait workflow
docs/                  model, inference, and reproducibility notes
```

## Reproducibility notes

- Use a new `OUT_ROOT` for every run. Task IDs begin at zero and existing output
  directories are not safe targets for a new manifest.
- Report rank, knot rule, knot seed, diagonal correction, likelihood convention,
  convergence, boundary estimates, and computing resources.
- Simulation seeds include the generating scenario and replicate, while all
  knot configurations within a replicate use the same simulated data.
- Do not interpret `c2` as a universal causal environmental effect. It is the
  variance associated with the supplied family/household relationship matrix.

## License and citation

The code is released under the BSD 3-Clause license. Add the repository URL and
the final manuscript citation to `CITATION.cff` after the GitHub repository and
article identifiers exist.
