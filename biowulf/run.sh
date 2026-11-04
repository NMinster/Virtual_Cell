#!/usr/bin/env bash
# =============================================================================
# Virtual_Cell :: run a stage in the RIGHT environment
# =============================================================================
# Thin wrapper that activates the `ovary_fm` conda env, sets the thread
# count from the SLURM allocation and exports the HPC proxy for the stages
# that talk to the network (fetch, fm, train, score --scorer cell-eval2).
# The pattern is Mouse_Compensatory/biowulf/run.sh.
#
# PREFER A BATCH JOB: bash biowulf/submit.sh <stage>. Run this directly
# only inside an sinteractive you are watching.
#
# USAGE (from the repo root):
#     bash biowulf/run.sh inspect
#     bash biowulf/run.sh fetch --dataset isola2024
#     bash biowulf/run.sh qc
#     bash biowulf/run.sh harmonise
#     bash biowulf/run.sh bench                       # the bundles and the ceilings table
#     bash biowulf/run.sh predict --model b1
#     bash biowulf/run.sh score --pred results/predictions/b1.h5ad --name b1
#     bash biowulf/run.sh t1
#     bash biowulf/run.sh fm --tool scprint2          # in the ovary_fm_scprint2 env, gpu partition
#     bash biowulf/run.sh train --arm c               # scVI/scANVI, in the ovary_fm_scvi env
#     bash biowulf/run.sh train --arm a               # writes the recipe only
#     bash biowulf/run.sh report
#     bash biowulf/run.sh synth --out /data/$USER/ovary_fm/fake
#     bash biowulf/run.sh python scripts/scbasecount_metadata.py --species mouse
#     bash biowulf/run.sh pip install -e .[census]    # install INTO the env
#     bash biowulf/run.sh pytest -q
#
# The roots come from the environment:
#     export OVARY_FM_DATA=/data/$USER/ovary_fm/data
#     export OVARY_FM_RESULTS=/data/$USER/ovary_fm/results
# =============================================================================
# CRLF guard: a checkout made on Windows (core.autocrlf=true) arrives with \r line
# endings and bash dies at `set -euo pipefail` with "set: pipefail\r: invalid option
# name". Re-run this script's text with the \r stripped; the permanent fix is
# `sed -i 's/\r$//' biowulf/*.sh` (the .gitattributes keeps fresh clones LF).
if [[ -z "${OVFM_CRLF_CLEANED:-}" ]] && grep -q $'\r' "$0" 2>/dev/null; then echo ">> $0 has CRLF line endings; running a cleaned copy (fix: sed -i 's/\r$//' biowulf/*.sh)" >&2; OVFM_CRLF_CLEANED=1 exec bash -c "$(tr -d '\r' < "$0")" "$0" "$@"; fi  # must stay one line ending in a comment
set -euo pipefail
: "${USER:=$(id -un)}"   # /data/$USER paths below; a bare shell may not export USER

ENV_NAME="${ENV_NAME:-ovary_fm}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

find_conda_base() {
    local c r n
    if [[ -n "${CONDA_INSTALL:-}" && -f "${CONDA_INSTALL}/etc/profile.d/conda.sh" ]]; then
        echo "${CONDA_INSTALL}"; return 0
    fi
    if command -v conda >/dev/null 2>&1; then
        c="$(conda info --base 2>/dev/null || true)"
        [[ -n "${c}" && -f "${c}/etc/profile.d/conda.sh" ]] && { echo "${c}"; return 0; }
    fi
    for c in "${CONDA_EXE:-}" "${MAMBA_EXE:-}"; do
        [[ -n "${c}" ]] || continue
        c="$(cd "$(dirname "${c}")/.." 2>/dev/null && pwd || true)"
        [[ -n "${c}" && -f "${c}/etc/profile.d/conda.sh" ]] && { echo "${c}"; return 0; }
    done
    for r in "${HOME}" "/data/${USER}" "/data/${USER}/conda"; do
        for n in miniforge3 mambaforge miniconda3 miniconda anaconda3 conda; do
            [[ -f "${r}/${n}/etc/profile.d/conda.sh" ]] && { echo "${r}/${n}"; return 0; }
        done
    done
    return 1
}

if base="$(find_conda_base)"; then
    # shellcheck disable=SC1091
    source "${base}/etc/profile.d/conda.sh"
else
    echo "!! No conda install found. Build the env first: bash biowulf/setup_env.sh" >&2
    exit 1
fi
set +u
conda activate "${ENV_NAME}" 2>/dev/null || {
    set -u
    echo "!! Could not activate '${ENV_NAME}'. Build it first: bash biowulf/setup_env.sh" >&2
    exit 1
}
set -u
if ! python -c "import scanpy" 2>/dev/null; then
    echo "!! '${ENV_NAME}' is active but scanpy is missing. Re-run: bash biowulf/setup_env.sh" >&2
    exit 1
fi
echo ">> env    : ${ENV_NAME}"
echo ">> python : $(command -v python)"

# weights and the Census cache off the home quota
export HF_HOME="${HF_HOME:-/data/${USER}/hf}"

# threads: match the numeric libraries to the CPU allocation
: "${OVARY_FM_THREADS:=${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 1)}}"
export OVARY_FM_THREADS
export OMP_NUM_THREADS="${OVARY_FM_THREADS}" OPENBLAS_NUM_THREADS="${OVARY_FM_THREADS}" \
       MKL_NUM_THREADS="${OVARY_FM_THREADS}" NUMEXPR_NUM_THREADS="${OVARY_FM_THREADS}"
echo ">> threads: ${OVARY_FM_THREADS}"

STAGE="${1:-}"
if [[ -z "${STAGE}" ]]; then
    echo "usage: bash biowulf/run.sh <stage|python|pip|pytest|synth> [args...]" >&2
    exit 1
fi
shift

# the proxy, for the stages that reach the network
case "${STAGE}" in
    fetch|fm|train|score|pip|python)
        export http_proxy="${http_proxy:-http://dtn03-e0:3128}"
        export https_proxy="${https_proxy:-http://dtn03-e0:3128}"
        ;;
esac

cd "${REPO_ROOT}"
case "${STAGE}" in
    python) exec python "$@" ;;
    pip)    exec python -m pip "$@" ;;
    pytest) exec python -m pytest "$@" ;;
    synth)  exec python scripts/make_synthetic.py "$@" ;;
    *)      exec ovary-fm "${STAGE}" "$@" ;;
esac
