#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# install_python_env.sh — create the analysis/plotting conda environment
# ─────────────────────────────────────────────────────────────────────────────
# Creates the conda env (groMD_env by default) from environment.yml. Run this
# once per cluster, on a login node — no GPU or SLURM job needed:
#
#   bash installation_scripts/install_python_env.sh
#
# The env name and conda module come from site_config.sh (GROMD_ENV, CONDA_MODULE),
# so this stays in sync with what the pipeline activates at run time.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../site_config.sh"

ENV_YML="${SCRIPT_DIR}/environment.yml"
[[ -f "$ENV_YML" ]] || { echo "[ERROR] environment.yml not found: $ENV_YML"; exit 1; }

module load "$CONDA_MODULE"

echo "[INFO] Creating conda env '${GROMD_ENV}' from ${ENV_YML} ..."
# -n overrides the name: field so the env always matches GROMD_ENV in site_config.
mamba env create -n "$GROMD_ENV" -f "$ENV_YML"

echo "[OK] Created conda env '${GROMD_ENV}'."
echo "     The pipeline activates it automatically (activate_python_env in site_config.sh)."
echo "     To use the analysis tools manually:  conda activate ${GROMD_ENV}"
echo ""
echo "     If the env already exists and you want to update it:"
echo "       conda env update -n ${GROMD_ENV} -f ${ENV_YML} --prune"
