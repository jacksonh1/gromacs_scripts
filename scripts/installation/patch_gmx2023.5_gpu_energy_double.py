#!/usr/bin/env python3
"""Make GROMACS 2023.5's CUDA nonbonded ENERGY accumulators double precision (forces untouched).

Usage:  patch_gmx2023.5_gpu_energy_double.py <gromacs-2023.5 source dir>

WHY: the CUDA nonbonded kernel reduces each warp's LJ / Coulomb energy and atomicAdd()s it
into ONE global `float` (eLJ, eElec). A float holding ~-3e5 kJ/mol resolves only ~0.02
kJ/mol, and thousands of warp additions land in a non-deterministic order, so the same
configuration evaluates to a total energy that differs by ~1 kJ/mol run to run. REST2's
PLUMED -hrex compares energies from separate evaluations, so every exchange carried
~0.37 kT of random noise (identical-lambda replicas accepted only ~80%; CPU is exact).
See knowledgebase/GOTCHAS.md, "on GPU, hrex exchange energies carry ~0.37 kT".

WHAT: the two global accumulators (device buffer + host staging copy) become `double`, and
the warp results are atomicAdd()ed as double. Per-pair energies and all FORCES are
unchanged, so trajectories are unchanged; only the reported nonbonded energy (and every
decision made from it) gets the precision back. The host still folds the value into
GROMACS' `real` (float) energy terms - the same precision a CPU run has.

REQUIRES compute capability >= 6.0 for atomicAdd(double*): build with
-DGMX_CUDA_TARGET_SM=89 (L40S). Only the CUDA backend is patched (the build is CUDA-only).

Every edit must match exactly once, or nothing is written.
"""
import sys
from pathlib import Path

EDITS = {
    # GROMACS builds its CUDA kernels with `-Xptxas -warn-double-usage -Xptxas -Werror`, a
    # developer guard against accidental double math. Our double atomicAdd is deliberate
    # (once per warp, energy steps only), so drop just that warning; -Werror stays for the rest.
    "cmake/gmxManageNvccConfig.cmake": [
        ("gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_FLAGS NVCC_HAS_PTXAS_WARN_DOUBLE_USAGE -Xptxas -warn-double-usage)\n",
         "# -warn-double-usage removed: nbnxm energy accumulators are double on purpose\n"
         "# (scripts/installation/patch_gmx2023.5_gpu_energy_double.py)\n"),
    ],
    "src/gromacs/nbnxm/gpu_types_common.h": [
        ("    //! LJ energy\n    float* eLJ = nullptr;\n    //! electrostatic energy\n    float* eElec = nullptr;",
         "    //! LJ energy (double: see patch_gmx2023.5_gpu_energy_double.py)\n    double* eLJ = nullptr;\n"
         "    //! electrostatic energy (double)\n    double* eElec = nullptr;"),
        ("    //! LJ energy output, size 1\n    DeviceBuffer<float> eLJ;\n    //! Electrostatics energy input, size 1\n    DeviceBuffer<float> eElec;",
         "    //! LJ energy output, size 1 (double: see patch_gmx2023.5_gpu_energy_double.py)\n    DeviceBuffer<double> eLJ;\n"
         "    //! Electrostatics energy input, size 1 (double)\n    DeviceBuffer<double> eElec;"),
    ],
    "src/gromacs/nbnxm/nbnxm_gpu_data_mgmt.cpp": [
        ("static_assert(sizeof(*nb->nbst.eLJ) == sizeof(float),",
         "static_assert(sizeof(*nb->nbst.eLJ) == sizeof(double),"),
        ("static_assert(sizeof(*nb->nbst.eElec) == sizeof(float),",
         "static_assert(sizeof(*nb->nbst.eElec) == sizeof(double),"),
    ],
    "src/gromacs/nbnxm/cuda/nbnxm_cuda_kernel.cuh": [
        ("    float*                     e_lj        = atdat.eLJ;\n    float*                     e_el        = atdat.eElec;",
         "    double*                    e_lj        = atdat.eLJ;\n    double*                    e_el        = atdat.eElec;"),
    ],
    "src/gromacs/nbnxm/cuda/nbnxm_cuda_kernel_utils.cuh": [
        ("reduce_energy_pow2(volatile float* buf, float* e_lj, float* e_el, unsigned int tidx)",
         "reduce_energy_pow2(volatile float* buf, double* e_lj, double* e_el, unsigned int tidx)"),
        ("        atomicAdd(e_lj, e1);\n        atomicAdd(e_el, e2);",
         "        atomicAdd(e_lj, static_cast<double>(e1));\n        atomicAdd(e_el, static_cast<double>(e2));"),
        ("reduce_energy_warp_shfl(float E_lj, float E_el, float* e_lj, float* e_el, int tidx, const unsigned int activemask)",
         "reduce_energy_warp_shfl(float E_lj, float E_el, double* e_lj, double* e_el, int tidx, const unsigned int activemask)"),
        ("        atomicAdd(e_lj, E_lj);\n        atomicAdd(e_el, E_el);",
         "        atomicAdd(e_lj, static_cast<double>(E_lj));\n        atomicAdd(e_el, static_cast<double>(E_el));"),
    ],
}
MARKER = ".pre_energy_double"


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    src = Path(sys.argv[1])
    new = {}
    for rel, edits in EDITS.items():
        path = src / rel
        text = path.read_text()
        if Path(str(path) + MARKER).exists():
            print(f"[SKIP] {rel} already patched ({MARKER} backup exists)")
            continue
        for old, rep in edits:
            n = text.count(old)
            if n != 1:
                raise SystemExit(f"[ERROR] {rel}: expected 1 match, found {n}:\n{old}")
            text = text.replace(old, rep)
        new[path] = text
    for path, text in new.items():
        Path(str(path) + MARKER).write_text(path.read_text())
        path.write_text(text)
        print(f"[OK] patched {path.relative_to(src)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
