# Per-target MOMA cohort view — design

**Date:** 2026-07-29
**Status:** approved, ready for implementation planning

## Problem

The pipeline renders one MOMA figure **per patient** (`viz/UMAP_MOMA_target_results_<id>.png`),
each showing that patient's top-10 gene knockouts. There is no way to ask the
complementary question: *for one gene, where does the whole cohort land when it is
knocked out?*

This design adds that inverted view — one figure **per target gene**, showing every
condition2 patient's wild-type position, their position after knocking out that gene,
and the nonIBD control cloud they are being compared against.

## Scope boundary (read first)

This is a **display** step. It computes no new statistic and makes no significance
claim. Per FUTURE_REQUESTS `[P0]`, the MOMA target ranking is hypothesis-only until
the competitive null is implemented, and the `[P0b]` boundary is explicit that
significance is never decided by 2D proximity in any embedding. The figure therefore
carries a caveat in its subtitle and must not be read as validating a target.

## Key enabling fact

No re-simulation is required. `moma_project.R` already projects **every gene** for
**every patient** into the predict-frame UMAP established by commit `e388856`, and
writes the result to `viz/UMAP_Results_MOMA_target_results_<id>.csv` with columns
`X`, `Y`, `wt_x`, `wt_y`, `recovery_proj_umap`, `recovery_proj_pca`. The control and
test baselines live in `viz/UMAP_predict_ref_df.csv`, in that same frame.

Re-slicing those files by gene is sufficient. Consequences:

- `moma_project.R` is **not modified** — the `[P0b]` fix stays untouched.
- No Gurobi, no MOMA re-run, no re-projection.
- The step works with no re-run on any already-completed run **whose MOMA
  projection was produced by the current predict-frame `moma_project.R`** — i.e.
  the two Human2 runs (`output_CD_ileum_Human2_tINIT_Human2_tINIT`,
  `output_UC_rectum_Human2_tINIT_Human2_tINIT`). The two Recon2M.2 runs
  (`output_CD_ileum_rbc_0.3`, `output_UC_rectum_rbc_0.3`) predate commit
  `e388856`: they have neither `viz/UMAP_predict_ref_df.csv` nor `wt_x`/`wt_y`
  columns, so the step **fails loud** there, as designed. Unlike `moma_viz.R`,
  which carries a legacy-frame fallback for those older results, this script
  deliberately has none — there is no correct frame to fall back to when the
  whole figure is a WT→KO displacement.
- Coordinates are by construction identical to the per-patient figures.

## Approach (chosen)

**Pure re-slice from existing outputs.** Two alternatives were rejected:

- *Re-project from `targeting_results/`* — reads ~5.8 GB across 36 patients to
  recompute numbers already on disk. Only needed for PC coordinates, which this
  figure does not use.
- *Extend `moma_project.R` to emit a tidy gene×patient table* — cheaper at plot time
  and a reusable artifact, but modifies a core file the `[P0b]` fix just stabilized,
  and every existing run would need the projection step re-run to gain the table.

## Components

All four changes are additive. No existing function's behaviour changes when the new
flag is absent.

| File | Change |
|---|---|
| `tf_metabolism/r_based_viz/moma_cohort_viz.R` | **New.** Renders one figure per target. |
| `tf_metabolism/utils.py` | **New** `visualize_cohort_ko_umap(output_dir, output_viz_dir, targets, mode='df')` — thin `subprocess` wrapper mirroring the R branch of `visualize_moma_umap`. |
| `tf_metabolism/utils.py` → `argument_parser` | **New flag** `--ko_targets`. |
| `tf_metabolism/__main__.py` | **One guarded call** + threading `ko_targets` through `run_pipeline`. |

### CLI flag

```
--ko_targets PKM ENO1 SLC33A1
--ko_targets PKM,ENO1,SLC33A1
```

`nargs='+'`, accepting space- **or** comma-separated values, matching the existing
`--model` / `--mem_method` convention. Default `None` → the step is skipped entirely,
so existing runs and reruns behave exactly as today.

### Gene resolution

A target may be given as an HGNC symbol (`PKM`) or a model gene ID
(`ENSG00000067225`). Resolution uses the run's own `viz/gene_map.csv` (columns
`entrez,symbol`; for Human2 runs the `entrez` column holds ENSG IDs — a known
misnomer, see CLAUDE.md §3h). Symbol matching is case-insensitive. This means the
same command works against a Recon2M.2 run and a Human2 run without the caller
knowing which namespace the run used.

If one symbol maps to several model IDs, every match is rendered, with the model ID
appended to the filename to keep outputs distinct.

**Resolution happens in the Python wrapper, not in R.** `visualize_cohort_ko_umap`
reads `gene_map.csv`, resolves each user-supplied name to one or more model IDs,
raises on an unresolvable name, and passes only **resolved model IDs** to the R
script. The R script still receives `--gene-map`, but uses it solely to label the
figure with the HGNC symbol — it performs no lookup of user input. This keeps the
one place that can reject user input in Python, where the error message can be
built with the same `pandas` machinery the rest of the pipeline uses, and leaves the
R script a pure renderer.

### R script interface

```
Rscript moma_cohort_viz.R \
    --output-dir     <run folder>        # for cohort_target_ranking.csv
    --output-viz-dir <run folder>/viz    # inputs and outputs
    --mode           df                  # default: df
    --targets        ENSG00000067225,ENSG00000074800   # resolved model IDs only
    --gene-map       <viz>/gene_map.csv                 # for figure labels only
```

Matches `moma_viz.R`'s existing `optparse` style and its `--gene-map` flag.

## Data flow

```
viz/UMAP_Results_MOMA_target_results_<id>.csv   (one per patient, gene-indexed)
        |  row for target gene -> X, Y, wt_x, wt_y, recovery_proj_umap
        v
viz/UMAP_predict_ref_df.csv                     (control + test baselines)
        |
        v
viz/ko_cohort_<SYMBOL>.png / .svg               (one per target)
```

## The figure

Per target gene, `viz/ko_cohort_<SYMBOL>.png` and `.svg`:

- **Control cloud** — nonIBD baselines from `UMAP_predict_ref_df.csv`, plus a cross
  at the control centroid. Colour `#009E73`, matching `moma_viz.R`.
- **Per patient** — an open WT dot at `(wt_x, wt_y)`, a filled KO dot at `(X, Y)`,
  and an arrow from WT to KO. One arrow per condition2 patient.
- **Arrow colour** — toward-control vs away-from-control, so the split is readable
  without counting arrows.
- **Subtitle** — `k/n patients moved toward control`, the count of patients dropped
  for lacking the gene, the gene's `rra_fdr` / `recurrence_prop` from
  `cohort_target_ranking.csv` when that file exists, and the `[P0]` hypothesis-only
  caveat.

### Toward vs away is decided by `recovery_proj_umap`, not `recovery_proj_pca`

The figure draws the UMAP frame, so the arrow colouring uses the signed projection
computed **in that same frame** (`recovery_proj_umap`). Using `recovery_proj_pca`
would colour arrows by a quantity the viewer cannot see, producing arrows that
visibly point toward control while coloured "away" (and vice versa).

`recovery_proj_pca` remains the scoring quantity everywhere else — this choice is
local to the display and does not touch `cohort_target_ranking`. It is a direct
application of the `[P0b]` display-vs-significance boundary.

## Error handling

Fail loud where the user explicitly asked for something that cannot be delivered;
warn and continue where part of the request is still satisfiable.

| Condition | Behaviour |
|---|---|
| `UMAP_predict_ref_df.csv` missing | **Fail loud** — the frame is undefined. |
| Zero `UMAP_Results_*.csv` found | **Fail loud** — MOMA projection has not run. |
| Gene name resolves to no model ID | **Fail loud**, listing near-miss candidates. |
| Target absent from *some* patients | Drop those patients; report the count in the subtitle and log. |
| Target absent from *all* patients | **Warn and skip** that gene; continue with the rest. |
| `cohort_target_ranking.csv` missing | Subtitle omits the ranking stats; figure still renders. |
| `Rscript` missing | **Fail loud** — consistent with every other MOMA step. |

The "absent from some patients" case is not an edge case: tINIT builds a
context-specific model per patient, so gene coverage genuinely varies. Measured on
the CD-ileum Human2 run: 2,041 genes common to all 36 patients, 2,285 in the union.

## Verification

The two completed **Human2** runs are the test fixture — verification needs no
re-run. (The two Recon2M.2 runs cannot serve as fixtures: they lack the
predict-frame reference, so the step correctly fails loud on them.)

1. `PKM` on CD-ileum Human2 (present in 36/36): 36 arrows, 36 WT dots, 36 KO dots.
2. WT dot coordinates match `wt_x`/`wt_y` in the per-patient CSVs, confirming the
   new figure shares the per-patient figures' frame.
3. `SLC33A1` (the one CD/UC shared target) renders in both cohorts.
4. A gene present in only a subset of patients: arrow count equals the present
   count, and the dropped count appears in the subtitle.
5. A bogus name (`NOTAGENE`) fails loud rather than producing an empty figure.
6. Flag omitted → no `ko_cohort_*` files, and the rest of the run is byte-identical
   to a pre-change run.

## Out of scope

- PCA-space variant of the figure (the UMAP frame was chosen deliberately).
- Auto-selecting targets from `cohort_target_ranking.csv` (the flag is opt-in only).
- Any change to how targets are scored or called significant — that is `[P0]`.
