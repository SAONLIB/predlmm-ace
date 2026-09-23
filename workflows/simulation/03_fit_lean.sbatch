#!/usr/bin/env bash
#SBATCH --job-name=sim_predlmm
#SBATCH --output=logs/03_fit_%A_%a.out
#SBATCH --error=logs/03_fit_%A_%a.err

set -euo pipefail
source "${WORKFLOW_DIR}/config.sh"
export PYTHONUNBUFFERED=1

PRIMARY_SE_METHOD="${PRIMARY_SE_METHOD:-structure}"
PRIMARY_SE_RANK="${PRIMARY_SE_RANK:-5000}"

row=$(awk -F'\t' -v id="${SLURM_ARRAY_TASK_ID}" \
    'NR > 1 && $1 == id {print; exit}' \
    "${OUT_ROOT}/manifests/factors.tsv")

if [[ -z "${row}" ]]; then
    echo "No factor manifest row for ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

IFS=$'\t' read -r fid rid sim n h c rep frac rank method <<< "${row}"
rdir=$(printf '%s/replicates/task_%06d' "${OUT_ROOT}" "${rid}")
fdir=$(printf '%s/factors/factor_%06d' "${OUT_ROOT}" "${fid}")
base=$(printf '%s/results/predlmm_%06d' "${OUT_ROOT}" "${fid}")
mkdir -p "${OUT_ROOT}/results"

continuous="PC1"
for ((j=2; j<=N_PCS; j++)); do
    continuous+=",PC${j}"
done

echo "factor=${fid} replicate_task=${rid} simulation=${sim} N=${n} h2=${h} c2=${c} replicate=${rep}"
echo "knot_method=${method} rank=${rank} primary_se=${PRIMARY_SE_METHOD}:${PRIMARY_SE_RANK}"

SECONDS=0
"${PYTHON_EXE}" -m predlmm_ace.fit \
    --factor-dir "${fdir}" \
    --metadata "${rdir}/metadata.tsv" \
    --phenotype phenotype \
    --family-column family_id \
    --continuous "${continuous}" \
    --model ACE \
    --out "${base}.json"
fit_seconds=${SECONDS}

"${PYTHON_EXE}" - "${base}.json" "${fit_seconds}" <<'PY'
import json, sys
path, seconds = sys.argv[1], float(sys.argv[2])
data = json.load(open(path))
data["fit_runtime_sec"] = seconds
data["runtime_sec"] = seconds
with open(path, "w") as handle:
    json.dump(data, handle, indent=2)
    handle.write("\n")
PY

if [[ "${method}" == "${PRIMARY_SE_METHOD}" && "${rank}" == "${PRIMARY_SE_RANK}" ]]; then
    echo "Calculating profile-Hessian SE for primary configuration"
    if "${PYTHON_EXE}" -m predlmm_ace.profile_se \
        --fit-json "${base}.json" \
        --factor-dir "${fdir}" \
        --metadata "${rdir}/metadata.tsv" \
        --phenotype phenotype \
        --family-column family_id \
        --continuous "${continuous}" \
        --out "${base}_se.json"; then
        total_seconds=${SECONDS}
        "${PYTHON_EXE}" - "${base}_se.json" "${fit_seconds}" "${total_seconds}" <<'PY'
import json, sys
path = sys.argv[1]
data = json.load(open(path))
data["fit_runtime_sec"] = float(sys.argv[2])
data["runtime_sec"] = float(sys.argv[3])
with open(path, "w") as handle:
    json.dump(data, handle, indent=2)
    handle.write("\n")
PY
    else
        echo "WARNING: point estimate completed but Hessian SE failed" >&2
    fi
else
    echo "Skipping SE for non-primary configuration"
fi

echo "COMPLETE PredLMM factor=${fid}"
