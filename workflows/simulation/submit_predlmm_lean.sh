#!/usr/bin/env bash
set -euo pipefail

WORKFLOW_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export WORKFLOW_DIR
requested_out_root="${OUT_ROOT:-}"
source "${WORKFLOW_DIR}/config.sh"

if [[ -n "${requested_out_root}" ]]; then
    export OUT_ROOT="${requested_out_root}"
else
    export OUT_ROOT="${WORKFLOW_DIR}/output_n12000_predlmm"
fi
export LEAN_REPS="${LEAN_REPS:-20}"
export LEAN_N="${LEAN_N:-12000}"
export FACTOR_SPECS="${FACTOR_SPECS:-structure:3000,structure:4000,structure:5000,random:5000}"
export PRIMARY_SE_METHOD="${PRIMARY_SE_METHOD:-structure}"
export PRIMARY_SE_RANK="${PRIMARY_SE_RANK:-5000}"
export LEAN_MAX_SIMULTANEOUS="${LEAN_MAX_SIMULTANEOUS:-20}"

dry_run=0
if [[ "${1:-}" == "--dry-run" ]]; then
    dry_run=1
elif [[ -n "${1:-}" ]]; then
    echo "Usage: bash submit_predlmm_lean.sh [--dry-run]" >&2
    exit 2
fi

for required in \
    config.sh make_lean_manifests.py \
    00_simulate.sbatch 01_make_grm.sbatch 01_knots.sbatch \
    02_factor.sbatch 03_fit_lean.sbatch 05_summary_lean.sbatch \
    summarize_lean.py; do
    test -s "${WORKFLOW_DIR}/${required}" || {
        echo "Missing ${WORKFLOW_DIR}/${required}" >&2
        exit 1
    }
done

"${PYTHON_EXE}" -c 'import numpy, pandas, scipy, predlmm_ace' || {
    echo "PYTHON_EXE lacks numpy/pandas/scipy or the installed predlmm_ace package: ${PYTHON_EXE}" >&2
    exit 1
}
test -x "${GCTA_EXE}" || command -v "${GCTA_EXE}" >/dev/null || {
    echo "Cannot execute GCTA_EXE=${GCTA_EXE}" >&2
    exit 1
}

mkdir -p "${OUT_ROOT}"/{logs,manifests,replicates,factors,results}
"${PYTHON_EXE}" "${WORKFLOW_DIR}/make_lean_manifests.py" \
    --reps "${LEAN_REPS}" \
    --N "${LEAN_N}" \
    --factor-specs "${FACTOR_SPECS}" \
    --out "${OUT_ROOT}/manifests"

nrep=$(( $(wc -l < "${OUT_ROOT}/manifests/replicates.tsv") - 1 ))
nfac=$(( $(wc -l < "${OUT_ROOT}/manifests/factors.tsv") - 1 ))

echo "Output:                  ${OUT_ROOT}"
echo "N:                       ${LEAN_N}"
echo "Replicates per scenario: ${LEAN_REPS}"
echo "Replicate tasks:         ${nrep}"
echo "Factor/fit tasks:        ${nfac}"
echo "Factor specifications:   ${FACTOR_SPECS}"
echo "SE configuration:        ${PRIMARY_SE_METHOD}:${PRIMARY_SE_RANK}"
echo "Concurrent array tasks:  ${LEAN_MAX_SIMULTANEOUS}"
echo "External methods:        NOT SUBMITTED"

if [[ "${FACTOR_SPECS}" == "structure:5000" && "${LEAN_REPS}" == "100" ]]; then
    [[ "${nrep}" -eq 600 && "${nfac}" -eq 600 ]] || {
        echo "Primary-run invariant failed: expected 600 replicate and factor tasks" >&2
        exit 1
    }
fi

if (( dry_run )); then
    echo "Dry run completed; no jobs submitted"
    exit 0
fi

opts=(
    --parsable
    --cpus-per-task="${SLURM_CPUS}"
    --mem="${SLURM_MEM}"
    --time="${SLURM_TIME}"
)
[[ -n "${SLURM_ACCOUNT}" ]] && opts+=(--account="${SLURM_ACCOUNT}")
[[ -n "${SLURM_PARTITION}" ]] && opts+=(--partition="${SLURM_PARTITION}")

submit_chunks() {
    local count=$1
    local script=$2
    local dependency=${3:-}
    local start=0 end array_spec job_id
    local ids=()

    while (( start < count )); do
        end=$((start + 899))
        (( end >= count )) && end=$((count - 1))
        array_spec="${start}-${end}%${LEAN_MAX_SIMULTANEOUS}"
        args=("${opts[@]}" --array="${array_spec}")
        [[ -n "${dependency}" ]] && args+=(--dependency="afterok:${dependency}")
        job_id=$(cd "${OUT_ROOT}" && sbatch "${args[@]}" "${script}")
        ids+=("${job_id}")
        start=$((end + 1))
    done

    local joined
    joined=$(IFS=:; echo "${ids[*]}")
    echo "${joined}"
}

simulate_jobs=$(submit_chunks "${nrep}" "${WORKFLOW_DIR}/00_simulate.sbatch")
grm_jobs=$(submit_chunks "${nrep}" "${WORKFLOW_DIR}/01_make_grm.sbatch" "${simulate_jobs}")
knot_jobs=$(submit_chunks "${nfac}" "${WORKFLOW_DIR}/01_knots.sbatch" "${simulate_jobs}")
factor_jobs=$(submit_chunks "${nfac}" "${WORKFLOW_DIR}/02_factor.sbatch" "${grm_jobs}:${knot_jobs}")
fit_jobs=$(submit_chunks "${nfac}" "${WORKFLOW_DIR}/03_fit_lean.sbatch" "${factor_jobs}")

summary_job=$(
    cd "${OUT_ROOT}"
    sbatch "${opts[@]}" --dependency="afterany:${fit_jobs}" \
        "${WORKFLOW_DIR}/05_summary_lean.sbatch"
)

printf 'stage\tjob_ids\n' > "${OUT_ROOT}/workflow_jobs.tsv"
printf 'simulate\t%s\n' "${simulate_jobs}" >> "${OUT_ROOT}/workflow_jobs.tsv"
printf 'grm\t%s\n' "${grm_jobs}" >> "${OUT_ROOT}/workflow_jobs.tsv"
printf 'knots\t%s\n' "${knot_jobs}" >> "${OUT_ROOT}/workflow_jobs.tsv"
printf 'factor\t%s\n' "${factor_jobs}" >> "${OUT_ROOT}/workflow_jobs.tsv"
printf 'predlmm\t%s\n' "${fit_jobs}" >> "${OUT_ROOT}/workflow_jobs.tsv"
printf 'summary\t%s\n' "${summary_job}" >> "${OUT_ROOT}/workflow_jobs.tsv"

cat "${OUT_ROOT}/workflow_jobs.tsv"
