#!/usr/bin/env bash
# =============================================================================
# Virtual_Cell :: one-shot Biowulf environment setup
# =============================================================================
# Creates (or updates) the `ovary_fm` conda env from environment.yml and
# installs the package editable. The other stacks are separate envs, on
# purpose (they pin different torch versions and would break each other's
# solve):
#
#     sinteractive --mem=16g --cpus-per-task=4     # ONLY if not already in a job
#     bash biowulf/setup_env.sh                    # main env + editable install
#     bash biowulf/setup_env.sh --arc              # + ovary_fm_arc (cell-eval2, Arc's scorer)
#     bash biowulf/setup_env.sh --train scprint2   # + ovary_fm_scprint2 (arm A and the fm stage)
#     bash biowulf/setup_env.sh --train tahoe      # + ovary_fm_tahoe (arm B)
#     bash biowulf/setup_env.sh --train scvi       # + ovary_fm_scvi (arm C, the baseline)
#     bash biowulf/setup_env.sh --fm uce           # + ovary_fm_uce (zero-shot arm)
#     bash biowulf/setup_env.sh --fm transcriptformer
#
# Idempotent. Overrides: ENV_NAME, CONDA_ROOT=/data/$USER/conda, CONDA_INSTALL.
# =============================================================================
# CRLF guard: a checkout made on Windows (core.autocrlf=true) arrives with \r line
# endings and bash dies at `set -euo pipefail` with "set: pipefail\r: invalid option
# name". Re-run this script's text with the \r stripped; the permanent fix is
# `sed -i 's/\r$//' biowulf/*.sh` (the .gitattributes keeps fresh clones LF).
if [[ -z "${OVFM_CRLF_CLEANED:-}" ]] && grep -q $'\r' "$0" 2>/dev/null; then echo ">> $0 has CRLF line endings; running a cleaned copy (fix: sed -i 's/\r$//' biowulf/*.sh)" >&2; OVFM_CRLF_CLEANED=1 exec bash -c "$(tr -d '\r' < "$0")" "$0" "$@"; fi  # must stay one line ending in a comment
set -euo pipefail
: "${USER:=$(id -un)}"   # /data/$USER paths below; a bare shell may not export USER

ENV_NAME="${ENV_NAME:-ovary_fm}"
CONDA_ROOT="${CONDA_ROOT:-/data/${USER}/conda}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
WITH_ARC=0; TRAIN=(); FM=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --arc)   WITH_ARC=1 ;;
        --train) shift; TRAIN+=("${1:?--train needs scprint2 | tahoe | scvi}") ;;
        --fm)    shift; FM+=("${1:?--fm needs uce | transcriptformer | scprint2}") ;;
        *) echo "!! unknown option $1" >&2; exit 1 ;;
    esac
    shift
done

export http_proxy="${http_proxy:-http://dtn03-e0:3128}"
export https_proxy="${https_proxy:-http://dtn03-e0:3128}"

find_conda_base() {
    local c r n
    if [[ -n "${CONDA_INSTALL:-}" && -f "${CONDA_INSTALL}/etc/profile.d/conda.sh" ]]; then
        echo "${CONDA_INSTALL}"; return 0
    fi
    if command -v conda >/dev/null 2>&1; then
        c="$(conda info --base 2>/dev/null || true)"
        [[ -n "${c}" && -f "${c}/etc/profile.d/conda.sh" ]] && { echo "${c}"; return 0; }
    fi
    for r in "${HOME}" "/data/${USER}" "/data/${USER}/conda"; do
        for n in miniforge3 mambaforge miniconda3 miniconda anaconda3 conda; do
            [[ -f "${r}/${n}/etc/profile.d/conda.sh" ]] && { echo "${r}/${n}"; return 0; }
        done
    done
    return 1
}
if ! base="$(find_conda_base)"; then
    echo ">> no conda found; installing Miniforge3 under ${CONDA_ROOT}/miniforge3"
    mkdir -p "${CONDA_ROOT}"
    curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
        -o "${CONDA_ROOT}/miniforge.sh"
    bash "${CONDA_ROOT}/miniforge.sh" -b -p "${CONDA_ROOT}/miniforge3"
    base="${CONDA_ROOT}/miniforge3"
fi
# shellcheck disable=SC1091
source "${base}/etc/profile.d/conda.sh"
SOLVER="conda"; command -v mamba >/dev/null 2>&1 && SOLVER="mamba"
echo ">> conda base: ${base} (${SOLVER})"
# envs and the package cache on /data, never the home quota
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-${CONDA_ROOT}/pkgs}"
export CONDA_ENVS_PATH="${CONDA_ENVS_PATH:-${CONDA_ROOT}/envs}"

build_env() {   # name yml
    local name="$1" yml="$2"
    if conda env list | awk '{print $1}' | grep -qx "${name}"; then
        echo ">> updating ${name} from ${yml}"
        ${SOLVER} env update -n "${name}" -f "${yml}" --prune
    else
        echo ">> creating ${name} from ${yml}"
        ${SOLVER} env create -n "${name}" -f "${yml}"
    fi
}

# Arm A needs the dev clones: pretraining lives in the scPRINT-2 repo, not
# in the PyPI package (their README). Editable installs into the env.
SRC_ROOT="${SRC_ROOT:-/data/${USER}/ovary_fm/src}"
scprint2_dev_clone() {
    mkdir -p "${SRC_ROOT}"
    local repo
    for repo in cantinilab/scPRINT-2 jkobject/scDataLoader cantinilab/GRnnData jkobject/benGRN; do
        local name="${repo#*/}"
        if [[ -d "${SRC_ROOT}/${name}/.git" ]]; then
            (cd "${SRC_ROOT}/${name}" && git pull --ff-only) || true
        else
            git clone --depth 1 "https://github.com/${repo}" "${SRC_ROOT}/${name}"
        fi
    done
    set +u; conda activate ovary_fm_scprint2; set -u
    python -m pip install -e "${SRC_ROOT}/scPRINT-2[dev]" -e "${SRC_ROOT}/scDataLoader[dev]" \
        -e "${SRC_ROOT}/GRnnData[dev]" -e "${SRC_ROOT}/benGRN[dev]"
    set +u; conda deactivate; set -u
    cat <<MSG
>> scPRINT-2 dev clones under ${SRC_ROOT}. Once, in the ovary_fm_scprint2 env:
>>     lamin init --storage /data/${USER}/ovary_fm/lamin --name ovary_fm --modules bionty
>>     scprint2 easy_setup
>>     hf download jkobject/scPRINT small-v2.ckpt --local-dir \${HF_HOME}/scprint
MSG
}

build_env "${ENV_NAME}" "${REPO_ROOT}/environment.yml"
set +u; conda activate "${ENV_NAME}"; set -u
python -m pip install -e "${REPO_ROOT}[census]"      # the Census API (Isola 2024) is an extra of the package
python -c "import ovary_fm, scanpy; print('ovary_fm', ovary_fm.__version__, 'scanpy', scanpy.__version__)"
set +u; conda deactivate; set -u

[[ "${WITH_ARC}" -eq 1 ]] && build_env ovary_fm_arc "${REPO_ROOT}/environment-arc.yml"
for t in "${TRAIN[@]:-}"; do
    [[ -z "${t}" ]] && continue
    case "${t}" in
        scprint2) build_env ovary_fm_scprint2 "${REPO_ROOT}/environment-train-scprint2.yml"
                  scprint2_dev_clone ;;
        tahoe)    build_env ovary_fm_tahoe "${REPO_ROOT}/environment-train-tahoe.yml" ;;
        scvi)     build_env ovary_fm_scvi "${REPO_ROOT}/environment-scvi.yml" ;;
        *) echo "!! --train ${t}: scprint2 | tahoe | scvi" >&2; exit 1 ;;
    esac
done
for m in "${FM[@]:-}"; do
    [[ -z "${m}" ]] && continue
    case "${m}" in
        scprint2) build_env ovary_fm_scprint2 "${REPO_ROOT}/environment-train-scprint2.yml" ;;
        uce)      build_env ovary_fm_uce "${REPO_ROOT}/environment-fm-uce.yml"
                  mkdir -p "${REPO_ROOT}/data/fm"
                  [[ -d "${REPO_ROOT}/data/fm/UCE" ]] || git clone --depth 1 https://github.com/snap-stanford/UCE "${REPO_ROOT}/data/fm/UCE" ;;
        transcriptformer) build_env ovary_fm_transcriptformer "${REPO_ROOT}/environment-fm-transcriptformer.yml" ;;
        *) echo "!! --fm ${m}: scprint2 | uce | transcriptformer" >&2; exit 1 ;;
    esac
done
echo ">> done. Next: export OVARY_FM_DATA and OVARY_FM_RESULTS, then bash biowulf/run.sh inspect"
