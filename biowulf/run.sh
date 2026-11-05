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
#     bash biowulf/run.sh pip install -e '.[census]'  # install INTO the env
#     bash biowulf/run.sh pytest -q
#     bash biowulf/run.sh convert huang2025           # the .rds under data/raw/huang2025/ -> mtx + meta.csv (loads R)
#     bash biowulf/run.sh rscript R/rds_to_mtx.R <in.rds> <outdir>/   # any R script in the module R
#
# `module load R` is refused on the login node, so `convert` and `rscript`
# belong in a batch job: bash biowulf/submit.sh convert huang2025.
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

if [[ -z "${OVARY_FM_DATA:-}" || -z "${OVARY_FM_RESULTS:-}" ]]; then
    echo ">> OVARY_FM_DATA / OVARY_FM_RESULTS are not set in this shell: data and results go under" >&2
    echo ">> ${REPO_ROOT}/{data,results}. For /data/${USER}, put the two exports in ~/.bashrc:" >&2
    echo ">>     export OVARY_FM_DATA=/data/${USER}/ovary_fm/data; export OVARY_FM_RESULTS=/data/${USER}/ovary_fm/results" >&2
fi

STAGE="${1:-}"
if [[ -z "${STAGE}" ]]; then
    echo "usage: bash biowulf/run.sh <stage|python|pip|pytest|synth|convert|rscript> [args...]" >&2
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

load_r() {
    # the module R carries Seurat and SingleCellExperiment; refused on the login node (use a job)
    if command -v module >/dev/null 2>&1 || [[ -f /usr/share/lmod/lmod/init/bash ]]; then
        # shellcheck disable=SC1091
        [[ -f /usr/share/lmod/lmod/init/bash ]] && source /usr/share/lmod/lmod/init/bash
        set +u; module load R >/dev/null; set -u
    fi
    command -v Rscript >/dev/null 2>&1 || { echo "!! Rscript is not available: run this inside a batch job or sinteractive (module load R is refused on the login node)" >&2; exit 1; }
    echo ">> R      : $(command -v Rscript) ($(Rscript --version 2>&1 | head -1))"
}

cd "${REPO_ROOT}"
case "${STAGE}" in
    python)  exec python "$@" ;;
    pip)     exec python -m pip "$@" ;;
    pytest)  exec python -m pytest "$@" ;;
    synth)   exec python scripts/make_synthetic.py "$@" ;;
    rscript) load_r; exec Rscript "$@" ;;
    convert)
        # bash biowulf/run.sh convert <dataset> [file.rds]: the .rds in that dataset's raw folder -> an mtx folder
        ds="${1:?convert needs a dataset name}"; shift
        dir="${OVARY_FM_DATA:-${REPO_ROOT}/data}/raw/${ds}"
        rds="${1:-$(ls -1 "${dir}"/*.rds 2>/dev/null | head -1)}"
        [[ -n "${rds}" && -f "${rds}" ]] || { echo "!! no .rds under ${dir}" >&2; exit 1; }
        load_r
        echo ">> convert: ${rds} -> ${dir}/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz,meta.csv}"
        exec Rscript R/rds_to_mtx.R "${rds}" "${dir}/" ;;
    *)       exec ovary-fm "${STAGE}" "$@" ;;
esac
