#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/config.sh"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
"$PYTHON_EXE" "$WORKFLOW_DIR/simulate_replicate.py" --simulation 1 --N 120 --M 300 --h2 .25 --c2 .10 --replicate 1 --seed 42 --n-pcs 3 --out "$tmp/data"
"$PYTHON_EXE" "$WORKFLOW_DIR/make_test_grm.py" --data "$tmp/data"
"$PYTHON_EXE" "$WORKFLOW_DIR/select_knots.py" --metadata "$tmp/data/metadata.tsv" --rank 30 --method structure --seed 42 --out "$tmp/knots.tsv"
"$PYTHON_EXE" -m predlmm_ace.grm_to_nystrom --grm-prefix "$tmp/data/grm/full" --knots "$tmp/knots.tsv" --out "$tmp/factor"
"$PYTHON_EXE" -m predlmm_ace.fit --factor-dir "$tmp/factor" --metadata "$tmp/data/metadata.tsv" --phenotype phenotype --family-column family_id --continuous PC1,PC2,PC3 --model ACE --out "$tmp/fit.json"
"$PYTHON_EXE" -m predlmm_ace.profile_se --fit-json "$tmp/fit.json" --factor-dir "$tmp/factor" --metadata "$tmp/data/metadata.tsv" --phenotype phenotype --family-column family_id --continuous PC1,PC2,PC3 --out "$tmp/fit_se.json"
test -s "$tmp/fit.json"; test -s "$tmp/fit_se.json"; echo 'SCALABLE SMOKE TEST PASSED'
