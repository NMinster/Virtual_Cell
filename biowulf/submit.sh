#!/usr/bin/env bash
# =============================================================================
# Virtual_Cell :: submit a stage as a BATCH job (the Mouse_Compensatory pattern)
# =============================================================================
#     bash biowulf/submit.sh fetch                   # downloads, with the proxy, off the login node
#     bash biowulf/submit.sh qc
#     bash biowulf/submit.sh chain qc harmonise bench   # one job each, chained with afterok
#     bash biowulf/submit.sh first                   # qc harmonise bench predict score t1 report (two jobs)
#     bash biowulf/submit.sh fm --tool scprint2      # gpu partition
#     bash biowulf/submit.sh train --arm c           # gpu partition (scVI)
#     bash biowulf/submit.sh python scripts/scbasecount_metadata.py --species mouse
#     bash biowulf/submit.sh --dry-run first         # print the sbatch lines only
#
# Everything after the stage is passed to run.sh. The roots come from the
# environment and are baked into the job script at submission time:
#     export OVARY_FM_DATA=/data/$USER/ovary_fm/data
#     export OVARY_FM_RESULTS=/data/$USER/ovary_fm/results
#
#     squeue -u $USER; tail -f logs/<stage>_<jobid>.out; scancel <jobid>
# Resources per stage can be overridden for one submission:
#     SUBMIT_MEM=200g SUBMIT_CPUS=16 SUBMIT_TIME=36:00:00 bash biowulf/submit.sh bench
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOG_DIR="${REPO_ROOT}/logs"

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" || "${1:-}" == "-n" ]]; then
    DRY_RUN=1; shift
fi
STAGE="${1:-}"
if [[ -z "${STAGE}" ]]; then
    echo "usage: bash biowulf/submit.sh [--dry-run] <stage> [run.sh args...]" >&2
    echo "stages: inspect fetch qc harmonise bench predict score t1 fm fm-compare train report" >&2
    echo "        chain <stage> [<stage>...]   one job each, chained with afterok" >&2
    echo "        first                        qc harmonise bench, then predict score t1 report" >&2
    echo "        python <script> [args]       any script on a compute node, with the proxy" >&2
    exit 1
fi
shift

# --- resources, by the heaviest stage in the list --------------------------
# A corpus of 1-2 M cells x ~30 k genes is tens of GB in memory; `bench`
# and `score` run Wilcoxon tests per bundle and are CPU-bound; `fm` and
# `train` on a real model want a GPU; `fetch` and `inspect` are light.
resources() {
    local heavy=0 gpu=0 light=0 other=0 stage
    for stage in $1; do
        case "${stage}" in
            harmonise|bench|score) heavy=1 ;;
            fm)    if [[ " $2 " == *" stub "* ]]; then light=1; else gpu=1; fi ;;
            train) if [[ " $2 " == *" c "* ]]; then gpu=1; else light=1; fi ;;
            inspect|fetch|report|python|predict) light=1 ;;
            *) other=1 ;;
        esac
    done
    if   [[ "${gpu}" -eq 1 ]];   then echo "8 64g 24:00:00 gpu --gres=gpu:a100:1"
    elif [[ "${heavy}" -eq 1 ]]; then echo "8 128g 12:00:00 norm"
    elif [[ "${light}" -eq 1 && "${other}" -eq 0 ]]; then echo "4 32g 06:00:00 norm"
    else                              echo "8 64g 08:00:00 norm"
    fi
}

submit() {
    local name="$1" stages="$2" dep="$3"; shift 3
    read -r cpus mem time partition gres <<< "$(resources "${stages}" "$*")"
    cpus="${SUBMIT_CPUS:-${cpus}}"; mem="${SUBMIT_MEM:-${mem}}"
    time="${SUBMIT_TIME:-${time}}"; partition="${SUBMIT_PARTITION:-${partition}}"
    local extra=()
    [[ -n "${gres:-}" ]] && extra+=("${gres}")
    if [[ -n "${dep}" ]]; then
        if [[ "${dep}" == after* ]]; then extra+=(--dependency="${dep}"); else extra+=(--dependency="afterok:${dep}"); fi
    fi
    mkdir -p "${LOG_DIR}"
    local exports="" var
    for var in OVARY_FM_DATA OVARY_FM_RESULTS HF_HOME; do
        if [[ -n "${!var:-}" ]]; then
            exports+="export ${var}=$(printf '%q' "${!var}")"$'\n'
        fi
    done
    local quoted=""
    if [[ $# -gt 0 ]]; then
        quoted="$(printf '%q ' "$@")"
    fi
    local script
    script="$(cat <<EOF
#!/bin/bash
set -uo pipefail
cd "${REPO_ROOT}"
${exports}echo ">> job \${SLURM_JOB_ID} on \$(hostname): stages [${stages}]"
echo ">> \$(date)"
echo ">> data=\${OVARY_FM_DATA:-<unset: config default>}  results=\${OVARY_FM_RESULTS:-<unset: config default>}"
status=0
for stage in ${stages}; do
    echo
    echo "################ \${stage} ################"
    bash biowulf/run.sh "\${stage}" ${quoted}|| { status=\$?; break; }
done
echo
echo ">> finished \$(date) with status \${status}"
exit \${status}
EOF
)"
    local cmd=(sbatch --job-name="ovfm_${name}" --cpus-per-task="${cpus}" --mem="${mem}" --time="${time}"
               --partition="${partition}" --output="${LOG_DIR}/${name}_%j.out" "${extra[@]}" --parsable)
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        { echo "${cmd[*]} <<'JOB'"; echo "${script}"; echo "JOB"; echo; } >&2
        echo "DRY_${name}"; return 0
    fi
    local id
    id="$(printf '%s\n' "${script}" | "${cmd[@]}")"
    echo "submitted ${name} as job ${id} (${cpus} CPU, ${mem}, ${time}, ${partition}${dep:+, after ${dep}})" >&2
    echo "${id}"
}

case "${STAGE}" in
    first)
        j1="$(submit object "qc harmonise bench" "" "$@")"
        j2="$(submit results "predict score t1 report" "${j1}" "$@")"
        echo "first results: ${j1} -> ${j2}; then: tail -f logs/results_${j2}.out"
        ;;
    chain)
        prev=""
        for stage in "$@"; do
            prev="$(submit "${stage}" "${stage}" "${prev}")"
        done
        ;;
    *)
        submit "${STAGE}" "${STAGE}" "" "$@" >/dev/null
        ;;
esac
