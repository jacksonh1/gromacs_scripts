#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# install_python_env.sh — create the analysis/plotting conda environment
# ─────────────────────────────────────────────────────────────────────────────
# Creates the conda env (groMD_env by default) from environment.yml. Run this
# once per cluster, on a login node — no GPU or SLURM job needed:
#
#   bash scripts/installation/install_python_env.sh
#
# The env name and conda module come from site_config.sh (GROMD_ENV, CONDA_MODULE),
# so this stays in sync with what the pipeline activates at run time.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../site_config.sh"

ENV_YML="${SCRIPT_DIR}/environment.yml"
[[ -f "$ENV_YML" ]] || { echo "[ERROR] environment.yml not found: $ENV_YML"; exit 1; }

module load "$CONDA_MODULE"

echo "[INFO] Creating conda env '${GROMD_ENV}' from ${ENV_YML} ..."
# -n overrides the name: field so the env always matches GROMD_ENV in site_config.
mamba env create -n "$GROMD_ENV" -f "$ENV_YML"


# Install the gromd_analysis package (editable) into the env. This is what puts
# the gromd-* console scripts on PATH — run_analysis.sh calls them by name, not by
# path. Editable, so an edit under gromd_analysis/ takes effect with no reinstall.
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
echo "[INFO] Installing gromd_analysis (editable) into '${GROMD_ENV}' ..."
set +u
eval "$(conda shell.bash hook)"
conda activate "$GROMD_ENV"
set -u
pip install --no-deps -e "$REPO_ROOT"

echo "[OK] Created conda env '${GROMD_ENV}' with gromd_analysis installed."
echo "     Entry points: gromd-layout gromd-plot-xvg gromd-plot-dssp gromd-cluster"
echo "                   gromd-acceptance gromd-chain-index"
echo "     The pipeline activates the env automatically (activate_python_env in site_config.sh)."
echo "     To use the analysis tools manually:  conda activate ${GROMD_ENV}"
echo ""
echo "     If the env already exists and you want to update it:"
echo "       conda env update -n ${GROMD_ENV} -f ${ENV_YML} --prune"
echo "       pip install --no-deps -e ${REPO_ROOT}"
echo ""
echo "     No-install fallback (e.g. a node without network): the package is a flat"
echo "     layout, so PYTHONPATH=${REPO_ROOT} makes 'import gromd_analysis' work."
