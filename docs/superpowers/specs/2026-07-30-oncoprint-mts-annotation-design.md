# Oncoprint mTS annotation + concordance scatter — design

**Date:** 2026-07-30
**Status:** approved, ready for implementation planning
**Kind:** STANDALONE PROTOTYPE — not wired into `run_pipeline`, no core file modified

## Problem

Joining the two Human2 runs' `cohort_target_ranking.csv` against
`viz/prototype_rmta_mTS.csv` (this session) shows the oncoprint's significant gene set
and the rMTA-style transformation score disagree, and disagree **asymmetrically**:

| cohort | significant genes | below own cohort's mTS background median | below 10th pctile | above 98th pctile |
|---|---|---|---|---|
| CD-ileum | 21 | **12** | **12** | 3 |
| UC-rectum | 40 | 2 | 1 | **30** |

The worst case is CD glycolysis: `ENO1/ENO2/ENO3` carry `rra_fdr = 2.0e-89` and
`recurrence_prop = 1.00` (top-10 lever in all 36 patients) yet sit at the **1.05th
percentile** of their cohort's all-gene mTS background. `GAPDH`/`GAPDHS` (1.8th),
`NME1`/`NME2` (1.5th), `CMPK1` (1.6th) and `BPGM`/`PGAM1-4` (6.5th) behave the same way.
Only `SLC7A1/2/3` (98.5-98.7th) and `PKM`/`PKLR` (95.6th) clear the 98th percentile, and
CD's strongest gene by mTS excess — `SLCO2A1`, +0.063, positive in 36/36 — is **absent
from the 21 significant genes entirely**.

Mechanism, and it is the known `[P0]` residual rather than a new defect:
`recovery_proj_pca` is a *signed* projection but is **not normalized by total flux
change**, so perturbation magnitude still leaks into the score. mTS (rMTA Eq. 8) divides
by total change, so it cannot be won that way. Spearman between the pipeline's
`median_score` and mTS excess is only **+0.380 (CD)** / **+0.251 (UC)** across ~2,200
genes — related but substantially different quantities. CD's contradicted genes are
exactly the high-throughput glycolytic/nucleotide hubs (`median_score` 3.5-21.2, PKM
21.2, ENO 13.1) against UC's 0.6-4.3.

Nothing on the current oncoprints reveals any of this. A reader sees ENO1/2/3 as the
strongest rows on the CD figure.

## Goal

Make the discordance **visible on the existing figure** without changing which genes it
shows, and quantify the overall metric relationship, so the user can decide whether mTS
is worth wiring into the pipeline.

## Non-goals (explicit)

- **Not** an mTS-*filtered* oncoprint. Filtering rows by mTS while the cells still encode
  per-patient `recovery_proj_pca` ranks produces a figure whose membership rule and body
  encode contradicting metrics, and would drop CD to 3 genes. Rejected this session.
- **Not** an mTS-*ranked* oncoprint (RRA over per-patient TS ranks). That is the coherent
  alternative and remains available — `transformation_scores()` already returns the full
  genes x patients matrix at `prototype_rmta_ts.py:105` and only collapses it to medians
  at line 131, so dumping it is a one-line change. Deferred pending this prototype.
- **Not** a significance test. mTS here is median-only; no bTS/wTS (those need new MTA
  MILP runs per gene per patient). `[P0]` stays open; nothing here validates a target.
- **No** change to any file under `tf_metabolism/`.

## Deliverables

Per cohort, written into that run's `viz/`:

1. `PROTO_target_oncoprint_mts_<COHORT>.png` / `.svg` — the existing oncoprint plus a
   far-right diverging mTS-excess bar.
2. `PROTO_mts_concordance_scatter_<COHORT>.png` / `.svg` — all-gene scatter of pipeline
   rank vs mTS excess.

One table at repo root, spanning both cohorts:

3. `PROTO_oncoprint_mts_concordance.csv` — 61 rows (21 CD + 40 UC), one per displayed
   gene.

All names are `PROTO_`-prefixed so they cannot collide with pipeline output, matching the
convention `prototype_ko_cohort_calibrated.R` established.

## Implementation approach: copy, not import

`utils.visualize_target_oncoprint` writes its figure and returns nothing, and the
annotation bar requires editing the fixed `add_axes` fractions from *inside* the layout
block (`utils.py:1288-1294`) — unreachable from a caller. The prototype therefore carries
its own copy of the render body (~300 lines).

This duplication will drift from `utils.py`. Accepted: the prototype's sole purpose is a
go/no-go decision, and a copy guarantees zero risk to the pipeline. Same trade-off
`prototype_ko_cohort_calibrated.R` made against `moma_cohort_viz.R`.

## Scope: which runs

The two **Human2** runs only:

- `output_CD_ileum_Human2_tINIT_Human2_tINIT`
- `output_UC_rectum_Human2_tINIT_Human2_tINIT`

The two Recon2M.2 runs have no `viz/prototype_rmta_mTS.csv`. A missing mTS table **fails
loud** with the exact command to generate it (`python prototype_rmta_ts.py ...`), matching
the repo's fail-loud convention — it must not silently render an unannotated figure.

## Data flow

1. Load `cohort_target_ranking.csv`, `viz/moma_gene_scores_full.csv`,
   `viz/prototype_rmta_mTS.csv`.
2. **Min-patient filter, required.** Restrict the mTS table to genes present in
   `>= 50%` of `n_patients.max()` before computing anything. Without it, genes present in
   only 1-2 patients' models top the ranking (CD: DCT/GNMT/SARDH; UC:
   ADSS2/TNRC6B/SLC17A7 — WORKLOG 2026-07-29c).
3. `background` = median `mTS_median` over the filtered genes. `excess = mTS_median −
   background`. `excess_pctile` = position of the gene's `mTS_median` within the sorted
   filtered background, as a percentage.
4. **Row selection is byte-identical to `visualize_target_oncoprint`:** all
   `significant == True`; cap at 40 by `rra_score`; if none significant, fall back to top
   20 by `rra_score`. The prototype must annotate exactly the rows production shows — no
   new gene filtering anywhere in this script.
5. Rank matrix, column ordering, GPR bracketing, organelle label colours, DE-flux `*`:
   all unchanged from the copied body.
6. A displayed gene absent from the *filtered* mTS table gets a **blank slot with a light
   grey tick, not a zero-length bar**, and is counted in the figure footnote. A missing
   measurement must not read as "sits at background".

## Figure 1: annotated oncoprint

**Bar.** Far right, new axis sharing the heatmap y-axis, diverging about zero with the
zero line drawn. `#0072B2` (blue) for `excess > 0`, `#D55E00` (vermillion) for
`excess < 0` — Okabe-Ito, colourblind-safe, and `#0072B2` already means "KO" on the §3i
figure. Panel header prints the cohort background value (e.g. `bg +0.256`) so the bar is
interpretable standalone. Symmetric x-limits at `max(|excess|)` rounded up, so the sign
is read from direction alone.

Placement puts the recurrence bar and the mTS bar adjacent, so a disagreement is one
horizontal scan per row — the reason this beat placing the bar next to the gene label.

**Layout fractions.** From the current values, `ax_main` narrows and a new axis is
inserted:

| axis | current | new |
|---|---|---|
| `ax_left` (labels, brackets) | `[0.015, 0.085, 0.150, 0.70]` | unchanged |
| `ax_main` (heatmap) | `[0.175, 0.085, 0.620, 0.70]` | `[0.175, 0.085, 0.545, 0.70]` |
| `ax_top` (patient bars) | `[0.175, 0.800, 0.620, 0.130]` | `[0.175, 0.800, 0.545, 0.130]` |
| `ax_right` (recurrence) | `[0.800, 0.085, 0.140, 0.70]` | `[0.725, 0.085, 0.105, 0.70]` |
| `ax_mts` (**new**) | — | `[0.845, 0.085, 0.135, 0.70]` |
| `ax_cbar` (rank colorbar) | `[0.800, 0.840, 0.140, 0.028]` | `[0.725, 0.840, 0.105, 0.028]` |

`fig_w` floor raised `8.0` -> `9.0` (`max(9.0, 0.24 * m + 5.8)`) to absorb the extra
column without compressing the heatmap.

**Footnote.** One line: mTS is a prototype implementation of rMTA Eq. 8 (median only, no
bTS/wTS), hypothesis-only, `[P0]` open; plus the count of displayed genes lacking an mTS
value.

## Figure 2: concordance scatter

One per cohort. Answers what the bar structurally cannot — whether the metrics disagree
only among the significant genes or across the whole model.

- x = rank by `rra_score` (1 = best pipeline target), log scale given ~2,200 genes.
- y = mTS excess, zero line drawn.
- All filtered genes as small grey points; `significant == True` genes overplotted in the
  same blue/vermillion sign encoding as the bar and labelled.
- **Labelling has no optional dependency.** `adjustText` is NOT installed in this
  environment (verified 2026-07-30; available: matplotlib 3.3.4, scipy 1.5.3, pandas,
  numpy). Use a deterministic greedy declutter: sort labels by y, and when a label would
  land within a fixed minimum vertical gap of the previous one, push it down by that gap.
  Stable and seed-free, consistent with the repo's determinism convention.
- `rho_rank` printed in the panel (see below), plus the background median.

**Two distinct Spearman statistics — do not conflate them.** Both were computed this
session and they differ:

| name | pair | CD | UC |
|---|---|---|---|
| `rho_rank` | rank by `rra_score` vs rank by `−excess`, all filtered genes | `+0.329` | `+0.302` |
| `rho_score` | `median_score` vs `excess`, all filtered genes | `+0.380` | `+0.251` |

The scatter's x-axis is the `rra_score` rank, so the **panel prints `rho_rank`**. Both go
in the stdout summary, each labelled.

## Table: `PROTO_oncoprint_mts_concordance.csv`

Columns: `cohort`, `Gene`, `Symbol`, `gpr_group`, `de_flux`, `rra_fdr`,
`recurrence_prop`, `median_score`, `mTS_median`, `background`, `excess`, `excess_pctile`,
`n_patients_mTS`. Sorted by `cohort` then `rra_score` ascending, so row order matches
figure row order.

## stdout summary

Per cohort: n significant, count below background, count below 10th percentile, count
above 98th percentile, both `rho_rank` and `rho_score` (each labelled), and the background
median.

## Acceptance criteria

No unit tests for a prototype. The port is correct iff the emitted CSV reproduces the
numbers computed ad hoc this session:

- **CD:** background `+0.2563`; 21 significant; **12** below background; **12** below
  10th pctile; **3** above 98th pctile; `rho_rank` `+0.329`; `rho_score` `+0.380`.
- **UC:** background `−0.0327`; 40 significant; **2** below background; **1** below 10th
  pctile; **30** above 98th pctile; `rho_rank` `+0.302`; `rho_score` `+0.251`.
- Spot values: CD `ENO1` excess `−0.03197` at pctile `1.052`; CD `SLC7A2` `+0.03154` at
  `98.72`; UC `ACOX1` `+0.09041` at `99.95`; UC `CPT2` `−0.009661` at `7.909`.

And visually: ENO1/2/3 point left (vermillion) on the CD oncoprint; the UC
peroxisomal/beta-oxidation arm points right (blue).

If the CSV disagrees with those numbers, the port is wrong — do not adjust the expected
numbers to match the code.

## Known caveat to carry forward

The GPR-collinear group in UC (`HACD1-4`/`HIGD2A` at `GPR_113`, `ELOVL3`/`ELOVL7` at
`GPR_224`, `ELOVL4`) has *negative* `median_score` (`−0.04` to `−0.42`) but 98.6-99th
percentile mTS. Because GPR-collinear genes have an identical knockout effect, neither
metric can separate their individual contributions. The figure will show them as blue
bars; that is not evidence for the individual gene, and the annotation does not resolve
the passenger-target concern recorded in CLAUDE.md §4.
