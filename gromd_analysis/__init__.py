"""Analysis tools for the GROMACS T-REMD / MD / REST2 pipelines.

The Python half of the analysis layer: trajectory-agnostic parsers, metrics and
plots that never invoke GROMACS. The GROMACS-driving steps (PBC correction,
stripping, alignment, gmx rms/gyrate/rmsf/dssp) stay as shell scripts in
`scripts/analysis/` — they are `gmx` CLI calls, and keeping them in shell keeps
every step echoed as a command you can re-run by hand.

Entry points installed by this package:

    gromd-layout        parse a job directory   (layout.py)
    gromd-plot-xvg      plot any GROMACS .xvg   (xvg.py)
    gromd-plot-dssp     secondary-structure map (dssp.py)
    gromd-cluster       conformational clusters (clustering.py)
    gromd-acceptance    REMD exchange rates     (remd_log.py)
    gromd-chain-index   per-chain .ndx groups   (chains.py)

Library use starts at JobDir, which resolves an OUTDIR into typed, existing paths:

    from gromd_analysis import JobDir
    job = JobDir.detect(Path("path/to/output_T-REMD/my-run"))
    job.mode, job.xtc, job.n_chains
"""

from gromd_analysis.chains import Chain, parse_molecules, protein_chains
from gromd_analysis.layout import JobDir, JobLayoutError
from gromd_analysis.remd_log import ExchangeStats, RemdLogError
from gromd_analysis.xvg import XvgData, parse_xvg

__all__ = [
    "Chain",
    "ExchangeStats",
    "JobDir",
    "JobLayoutError",
    "RemdLogError",
    "XvgData",
    "parse_molecules",
    "parse_xvg",
    "protein_chains",
]
