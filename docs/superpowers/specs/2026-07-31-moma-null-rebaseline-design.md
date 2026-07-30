# MOMA null re-baselining — design

**Date:** 2026-07-31
**Status:** approved, implementing
**Kind:** STANDALONE PROTOTYPE — not wired into `run_pipeline`, no core file modified

## The question

Is `[P0]`'s "global drift" simply the projection of an unreachable wild-type anchor onto
MOMA's feasible set — an artifact shared by every knockout by construction — rather than a
property of the cohort?

## What is already measured (2026-07-31, this session)

Splitting each knockout into the part shared by all genes and the part specific to the gene,
for one representative patient per cohort:

| | shared by every gene | least-perturbing gene | gene-specific part | ratio |
|---|---|---|---|---|
| CD (H4031) | 2846 | 2842 | 47.6 | **59.7x** |
| UC (H4035) | 4232 | 4230 | 76.5 | **55.3x** |

- **Not a reference mismatch.** `flux2.csv` is bit-identical to MOMA's actual anchor
  (`flux2/Flux_prediction_<sample>.csv`) — max difference 0 across all 7,975 reactions,
  identical reaction sets. So the offset is not a bookkeeping error.
- **The mechanism is the anchor itself.** The reactions carrying the offset sit at the flux
  bounds — MAR09005 `-1000 -> -514.68`, MAR00477 `500 -> 63.95`, MAR04896 `-1000 -> -577.25`,
  MAR08104 `-1000 -> -592.24`, MAR02025 `1000 -> 593.60`. The flux2 prediction leaves many
  reactions pinned at +-1000; MOMA's quadratic objective penalises large magnitudes and
  shrinks them. The pull is broad, not a few outliers: the top 5 reactions are only 11.6%
  of the squared offset.
- **No gene reproduces the anchor.** 0 of 2,172 genes return the WT flux exactly, and the
  least-perturbing gene is 2842 away while the consensus offset is 2846.

**Why this matters beyond the figures.** mTS is computed on `delta = KO - WT`, so it contains
the offset. CD's background of `+0.2563` with 100% of genes positive may therefore be the
direction of this projection, shared by construction. That is an inference, not yet shown —
this prototype tests it.

## Design

**`prototype_moma_null_rebaseline.py`** (repo root, standalone, `PROTO_` outputs).

### One scoring function, two baselines
A single vectorised rMTA Eq. 8 implementation taking the baseline as a parameter, reusing
`prototype_rmta_ts.desired_direction()` for the per-reaction desired directions. Both
branches MUST be scored by identical code — scoring them separately would leave any
difference attributable to the implementation rather than the baseline.

### Stage A — re-baselined mTS (no solver, minutes)
Per patient: `null_r = median over all ~2,170 genes` of the KO flux for reaction `r`.
Compute `TS_wt` on `delta = KO - WT` and `TS_null` on `delta = KO - null`. Aggregate per gene
to medians across patients under the >=50%-of-patients filter (REQUIRED — genes present in
1-2 patients otherwise top the ranking, WORKLOG 2026-07-29c).

### Stage B — validate the consensus proxy (Gurobi, ~10 min)
The consensus is a proxy: every column it averages does have some gene knocked out. For the
first 5 patients per cohort by sample ID (deterministic), load the reconstructed SBML and run
`Simulator.run_MOMA(wild_flux=WT, flux_constraints={})` — nothing knocked out — to obtain the
true null. Report `||true_null - consensus||` against the gene-specific scale (47.6 CD /
76.5 UC) and the offset (2846 / 4232). Gurobi is licensed and all 36 CD / 20 UC reconstructed
models are on disk (verified 2026-07-31).

## Acceptance criteria

1. **Regression gate — non-negotiable.** With `baseline=WT` the script must reproduce the
   existing mTS numbers exactly: CD background `+0.2563`, 100% of 2,187 filtered genes
   positive; UC background `-0.0327`, 1.9% of 2,162 positive. If it does not, the
   reimplementation is unfaithful and nothing downstream is interpretable. Do NOT adjust the
   expected numbers to match the code.
2. **The test.** Does the CD background collapse from `+0.2563` toward 0, and does "100% of
   genes positive" fall toward ~50%?
3. **Does the ranking survive?** Spearman(old mTS excess, new `TS_null`) per cohort, plus the
   explicit ordering of SLC7A2 / ENO1 / ACOX1 / CPT2. This matters as much as the background:
   if re-baselining reshuffles which genes win, then every target list produced so far was
   scoring the projection rather than the knockout.
4. **Stage B.** Pass if `||true_null - consensus||` is small relative to the gene-specific
   scale; report the ratio either way rather than asserting a threshold post hoc.

## Outputs

- `PROTO_moma_null_rebaseline_<cohort>.csv` — per gene: `mTS_wt`, `mTS_null`, the percentile
  of each within its own background, `n_patients`.
- `PROTO_moma_null_validation.csv` — per patient: `||consensus - WT||`,
  `||true_null - consensus||`, `||true_null - WT||`, and the ratio to the gene-specific scale.
- One figure per cohort: the all-gene background distribution before vs after re-baselining,
  with SLC7A2 / ENO1 / ACOX1 / CPT2 marked.
- stdout: the two backgrounds, % positive, the Spearman, and the four genes' ordering.

## Non-goals

- Not wired into `run_pipeline`; changes no score the pipeline currently uses.
- Not a significance test, and still mTS-only (no bTS/wTS), so `[P0]` is not closed by a
  passing result — it would be *explained*, which is a different and prior thing.
- Does not attempt to fix the flux2 prediction. Whether the +-1000 pinning reflects genuine
  unbounded loops in the reconstructed models or a bounds-setting choice in the flux2 step is
  a separate question, unresolved, and noted here so it is not forgotten.

## Interpretation guide, written before the result

- **If the background collapses AND the ranking is preserved** — the drift is the projection;
  existing target lists survive but every absolute mTS number quoted so far was inflated by a
  constant, and the calibration story simplifies.
- **If the background collapses AND the ranking changes** — the more serious outcome. Every
  target list to date, including the oncoprint sets, was partly scoring the projection.
- **If the background does NOT collapse** — the drift is not the projection, `[P0]` stands as
  written, and the offset is a separate (still real) figure-level artifact.
