#!/usr/bin/env bash
#SBATCH --job-name=sim_predlmm
#SBATCH --output=logs/03_fit_%A_%a.out
#SBATCH --error=logs/03_fit_%A_%a.err
set -euo pipefail; source "${WORKFLOW_DIR}/config.sh"; export PYTHONUNBUFFERED=1
row=$(awk -F'\t' -v id="$SLURM_ARRAY_TASK_ID" 'NR>1&&$1==id{print;exit}' "$OUT_ROOT/manifests/factors.tsv"); IFS=$'\t' read -r fid rid sim n h c rep frac rank method <<<"$row"
rdir=$(printf '%s/replicates/task_%06d' "$OUT_ROOT" "$rid"); fdir=$(printf '%s/factors/factor_%06d' "$OUT_ROOT" "$fid"); base=$(printf '%s/results/predlmm_%06d' "$OUT_ROOT" "$fid"); mkdir -p "$OUT_ROOT/results"
cont=$(printf 'PC1'); for ((j=2;j<=N_PCS;j++)); do cont+=",PC${j}"; done
SECONDS=0; "$PYTHON_EXE" -m predlmm_ace.fit --factor-dir "$fdir" --metadata "$rdir/metadata.tsv" --phenotype phenotype --family-column family_id --continuous "$cont" --model ACE --out "${base}.json"
if "$PYTHON_EXE" -m predlmm_ace.profile_se --fit-json "${base}.json" --factor-dir "$fdir" --metadata "$rdir/metadata.tsv" --phenotype phenotype --family-column family_id --continuous "$cont" --out "${base}_se.json"; then
  "$PYTHON_EXE" -c "import json;p='${base}_se.json';z=json.load(open(p));z['runtime_sec']=$SECONDS;open(p,'w').write(json.dumps(z,indent=2)+'\\n')"
else
  echo 'WARNING: point estimate completed but ordinary Hessian SE was invalid; retaining point-estimate JSON' >&2
fi
