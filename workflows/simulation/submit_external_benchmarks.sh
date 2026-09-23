#!/usr/bin/env bash
set -euo pipefail

WORKFLOW_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export WORKFLOW_DIR
requested_out_root="${OUT_ROOT:-}"
source "${WORKFLOW_DIR}/config.sh"

if [[ -z "${requested_out_root}" ]]; then
    echo "Set OUT_ROOT to the completed lean simulation directory" >&2
    exit 1
fi
export OUT_ROOT="${requested_out_root}"

export RUN_GCTA="${RUN_GCTA:-1}"
export RUN_FASTGWA="${RUN_FASTGWA:-0}"
export GCTA_REPS="${GCTA_REPS:-20}"
export FASTGWA_REPS="${FASTGWA_REPS:-10}"
export FASTGWA_CUTOFF="${FASTGWA_CUTOFF:-0.05}"
export EXTERNAL_MAX_SIMULTANEOUS="${EXTERNAL_MAX_SIMULTANEOUS:-10}"

manifest="${OUT_ROOT}/manifests/replicates.tsv"
test -s "${manifest}" || { echo "Missing ${manifest}" >&2; exit 1; }

missing_grm=0
while IFS=$'\t' read -r id sim n h c rep; do
    [[ "${id}" == "task_id" ]] && continue
    prefix=$(printf '%s/replicates/task_%06d/grm/full' "${OUT_ROOT}" "${id}")
    [[ -s "${prefix}.grm.id" ]] || { echo "Missing ${prefix}.grm.id" >&2; missing_grm=1; }
done < "${manifest}"
(( missing_grm == 0 )) || exit 1

opts=(--parsable --cpus-per-task="${SLURM_CPUS}" --mem="${SLURM_MEM}" --time="${SLURM_TIME}")
[[ -n "${SLURM_ACCOUNT}" ]] && opts+=(--account="${SLURM_ACCOUNT}")
[[ -n "${SLURM_PARTITION}" ]] && opts+=(--partition="${SLURM_PARTITION}")

mkdir -p "${OUT_ROOT}/external" "${OUT_ROOT}/logs"
job_file="${OUT_ROOT}/external_jobs.tsv"
printf 'stage\tjob_id\n' > "${job_file}"

if [[ "${RUN_GCTA}" == "1" ]]; then
    gcta_ids=$(awk -F'\t' -v maxrep="${GCTA_REPS}" 'NR>1 && $6<=maxrep {ids=(ids?ids ",":"") $1} END{print ids}' "${manifest}")
    [[ -n "${gcta_ids}" ]] || { echo "No GCTA tasks selected" >&2; exit 1; }
    gcta_job=$(cd "${OUT_ROOT}" && sbatch "${opts[@]}" \
        --array="${gcta_ids}%${EXTERNAL_MAX_SIMULTANEOUS}" \
        "${WORKFLOW_DIR}/04_gcta_only.sbatch")
    printf 'gcta_unrelated\t%s\n' "${gcta_job}" >> "${job_file}"
fi

if [[ "${RUN_FASTGWA}" == "1" ]]; then
    fastgwa_ids=$(awk -F'\t' -v maxrep="${FASTGWA_REPS}" 'NR>1 && $6<=maxrep {ids=(ids?ids ",":"") $1} END{print ids}' "${manifest}")
    [[ -n "${fastgwa_ids}" ]] || { echo "No fastGWA tasks selected" >&2; exit 1; }
    fastgwa_job=$(cd "${OUT_ROOT}" && sbatch "${opts[@]}" \
        --array="${fastgwa_ids}%${EXTERNAL_MAX_SIMULTANEOUS}" \
        "${WORKFLOW_DIR}/04_fastgwa_only.sbatch")
    printf 'fastgwa_cutoff_%s\t%s\n' "${FASTGWA_CUTOFF}" "${fastgwa_job}" >> "${job_file}"
fi

cat "${job_file}"
echo "After external jobs finish, run: bash resummarize_lean.sh"
