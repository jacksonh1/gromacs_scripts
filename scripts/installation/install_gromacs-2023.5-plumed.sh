#!/bin/bash
#SBATCH -J build_gmx2023.5-plumed
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p mit_normal_gpu           # build on any FREE node (your pi_keating queue is busy)
#SBATCH --constraint="rocky8"       # binary must match the run nodes' OS (see CLAUDE.md gotcha)
#SBATCH --gres=gpu:l40s:1           # require L40S = same GPU as the run target (sm_89); used by make check
#SBATCH --mem=80000
#SBATCH -t 3:00:00
#SBATCH -o ./logs/build_gmx2023.5-plumed_%j.out
#SBATCH -e ./logs/build_gmx2023.5-plumed_%j.err
# ─────────────────────────────────────────────────────────────────────────────
# Build GROMACS 2023.5 patched with PLUMED 2.9.4, for REST2 / -hrex.
#
# WHY 2023.5 (not 2024.3): PLUMED's hrex energy machinery is NOT ported to the
# GROMACS-2024 patch — REST2 runs but silently gives ZERO exchanges. PLUMED dev
# G. Bussi: "hrex is not yet supported in the 2024 patch." 2023.5 is the newest
# GROMACS with working PLUMED hrex. See knowledgebase/plans/REST2-pipeline.md and
# the "-hrex silently broken" gotcha in CLAUDE.md. Your existing 2024.3-plumed
# build is UNTOUCHED and keeps serving T-REMD / MD.
#
# SIMD PORTABILITY: GROMACS bakes the build node's CPU SIMD/-march into the binary, so a
# binary built on a newer microarch can crash with "illegal instruction" on older run nodes.
# Your pi_keating run nodes (node[3619-3620]) are Skylake-AVX512. Rather than pin the build
# to those (busy) nodes, we build on ANY free node and FORCE Skylake-AVX512 via
# -DGMX_SIMD=AVX_512 + -march=skylake-avx512 (cmake call below). gcc emits Skylake-AVX512
# code regardless of the build host's CPU, so the binary runs on BOTH your pi_keating nodes
# and the newer mit_normal_gpu nodes. We still require an L40S GPU (same as the run target,
# sm_89) so make check + the CUDA arch selection match the real GPU. See the SIMD gotcha in CLAUDE.md.
#
# Mirrors install_gromacs-plumed.sh (same modules, plumed patch step); GMX_VER, install
# prefix, CUDA version, SIMD pinning, and build partition differ. Uses a pre-downloaded
# source tree at $SRC (no download step).
#
# Submit from a dir that has a ./logs/ subdir (e.g. the repo root):
#   sbatch scripts/installation/install_gromacs-2023.5-plumed.sh
# ─────────────────────────────────────────────────────────────────────────────
set -uo pipefail

GMX_VER=2023.5
PREFIX="$HOME/opt/gromacs/${GMX_VER}-plumed"          # parallel to 2024.3-plumed
SRC="${SRC:-/orcd/pool/004/jhalpin/installations/gromacs-2023.5}"   # pre-downloaded, clean, unpatched

# ── Environment ──────────────────────────────────────────────────────────────
# Same as the 2024.3 build EXCEPT CUDA: GROMACS 2023.5 predates cuda 12.9, so use
# cuda/12.4.0 (contemporaneous, well-supported) via the deprecated-modules tree.
# NOTE: the REST2 run engine must load the SAME cuda to match this binary
# (module load deprecated-modules; module load cuda/12.4.0).
export LD_LIBRARY_PATH=
module load openmpi/5.0.8
module load deprecated-modules
module load cuda/12.4.0
module load cmake/3.24.3-x86_64                        # default cmake is 4.x, which rejects GROMACS 2023.5's
                                                       # bundled googletest (cmake_minimum_required < 3.5); use era cmake
source ~/plumed.sh                                     # sets PLUMED_KERNEL / LD_LIBRARY_PATH (PLUMED 2.9.4)

echo "[INFO] host=$(hostname)"
echo "[INFO] plumed: $(command -v plumed)  ($(plumed info --version 2>/dev/null))"
echo "[INFO] target: $PREFIX"

# ── CUDA note ────────────────────────────────────────────────────────────────
# cuda/12.4.0 (loaded above via deprecated-modules) is a good match for GROMACS
# 2023.5, so no version-check trouble is expected. If for any reason it's needed,
# a CPU-only build (-DGMX_GPU=OFF) also works — hrex is correct on CPU (verified),
# just slower. Do NOT use cuda/12.9 or 13.x with 2023.5.

# ── Use the pre-downloaded source (no download) ──────────────────────────────
mkdir -p ./logs
[[ -d "$SRC" ]] || { echo "[ERROR] GROMACS source not found: $SRC"; exit 1; }
cd "$SRC" || exit 1

# ── PLUMED patch (must be done on the LOGIN node beforehand) ──────────────────
# The standalone `plumed` CLI does NOT run on the compute nodes here (missing runtime
# dependency — same class as the hcoll issue; the plumed *kernel* loads fine at mdrun
# time, just not the CLI). So `plumed patch` fails inside this batch job. The runtime
# ("shared") patch is only needed at PATCH time, not at build time, so we patch the
# source on the login node first and just build here. On the LOGIN node, once:
#     source ~/plumed.sh && cd "$SRC" && rm -rf build && plumed patch -p -e gromacs-2023.5
if ls "$SRC"/src/gromacs/mdrun/replicaexchange.cpp.preplumed >/dev/null 2>&1; then
  echo "[INFO] source is PLUMED-patched (.preplumed backups present) — skipping CLI patch"
else
  echo "[ERROR] source is NOT PLUMED-patched. On the LOGIN node run:"
  echo "          source ~/plumed.sh && cd $SRC && rm -rf build && plumed patch -p -e gromacs-${GMX_VER}"
  echo "        then resubmit this build."
  exit 1
fi
rm -rf build && mkdir build && cd build || exit 1

# ── Configure (2024.3 flags + forced Skylake-AVX512 SIMD/-march for portability) ─
cmake .. \
  -DGMX_THREAD_MPI=OFF \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DREGRESSIONTEST_DOWNLOAD=ON \
  -DGMX_GPU=CUDA \
  -DGMX_MPI=on \
  -DGMX_SIMD=AVX_512 \
  -DCMAKE_C_FLAGS="-march=skylake-avx512" \
  -DCMAKE_CXX_FLAGS="-march=skylake-avx512" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  || { echo "[ERROR] cmake failed — if it's a CUDA error, try the CPU-only fallback -DGMX_GPU=OFF (see CUDA note above)"; exit 1; }

# ── Build ────────────────────────────────────────────────────────────────────
make -j 8 || { echo "[ERROR] make failed"; exit 1; }

# make check is informative but can fail on unrelated GPU/precision regression
# tests; do NOT let that block install. The authoritative REST2 test is the
# accept8 hrex rerun below, not the regression suite.
make check || echo "[WARN] 'make check' reported failures — review, but proceeding to install (real validation is the hrex acceptance test)"

make install || { echo "[ERROR] make install failed"; exit 1; }

echo ""
echo "=========================================================================="
echo "[OK] Installed GROMACS ${GMX_VER}+PLUMED at: $PREFIX"
echo ""
echo "NEXT — verify hrex actually exchanges (should be NONZERO, unlike 2024.3):"
echo "  1) Rebuild the 8-replica tprs against this GROMACS, then rerun the"
echo "     staged fine-ladder acceptance test in"
echo "       /orcd/data/keating/001/jhalpin/MD/rest2_smoke/"
echo "     pointing GMXRC at:  ${PREFIX}/bin/GMXRC"
echo "  2) Expect ~20-40% acceptance (and P=1.0 for a scale-1.0 sanity pair)."
echo ""
echo "Do NOT change site_config.sh GMXRC globally — that would switch T-REMD/MD"
echo "onto 2023.5 too. The REST2 engine will point at this build explicitly."
echo "=========================================================================="
