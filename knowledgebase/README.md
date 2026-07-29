# gromacs_REMD Knowledgebase

A living record of the T-REMD / MD / REST2 GROMACS pipelines: what we've done,
why we decided what we decided, and the science behind the protocols. Maintained
by Claude across sessions.

## Files

| File | What it holds |
|------|---------------|
| [`HISTORY.md`](HISTORY.md) | Dated log — one entry per working session/day. Newest on top. |
| [`DECISIONS.md`](DECISIONS.md) | The *why* behind design choices. One entry per decision; link the code. |
| [`SCIENCE.md`](SCIENCE.md) | Protocol rationale — what each engine does, ensembles, when to use which. |
| [`plans/`](plans/) | Design/plan docs for in-progress features (e.g. `plans/REST2-pipeline.md`). |

## What lives elsewhere (not here — don't duplicate)

- **`CLAUDE.md`** (repo root) — the operational contract: working rules + the
  **"Known gotchas"** pitfalls database. Auto-loaded every session, which is what
  keeps the same mistake from happening twice, so gotchas **stay there**. This KB
  references them; it does not copy them.
- **`scripts/simulation/*-output-guide.md`**, **`scripts/analysis/*_reference.md`**,
  **`scripts/analysis/README.md`** — user-facing, per-artifact docs; they live next
  to the code they describe.
- **`README.md`** (repo root) — the user-facing project overview + feature list.

## Maintenance protocol (for Claude)

At the **end of a working session** where something non-trivial happened:

1. **`HISTORY.md`** — add a dated entry: what we worked on, what changed, what's
   still open. Terse and factual — a log, not a report.
2. **`DECISIONS.md`** — if a non-obvious choice was made, record the *why* (and the
   alternative rejected). Never fabricate rationale — if the *why* isn't known, say so.
3. **`SCIENCE.md`** — update if the scientific protocol or its rationale changed.
4. **A gotcha discovered?** → it goes in **`CLAUDE.md`**, not here.
5. **A feature graduating from a plan?** → move its design notes out of `plans/`
   into DECISIONS/SCIENCE + the user-facing guide, and leave a HISTORY entry.

Rules: **be terse** (low cognitive load is a project value). Prefer appending over
rewriting history. For "where is X", use `rg -n` / codegraph against the live tree
rather than a hand-maintained index.
