# Per-target MOMA Cohort View Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in pipeline step that renders one figure per knocked-out target gene, showing every condition2 patient's wild-type → knockout shift against the nonIBD control cloud.

**Architecture:** Pure re-slice of data already on disk. `moma_project.R` writes one CSV per patient containing every gene's projected knockout coordinates; the new step pivots those files by gene instead of by patient. A new R script renders, a thin Python wrapper resolves gene names and shells out to it, and `run_pipeline` calls the wrapper when `--ko_targets` is given. No MOMA re-run, no change to `moma_project.R`.

**Tech Stack:** Python 3.6.15 + pandas 0.23.4, R 4.2.0 + ggplot2 3.5.1 + optparse 1.7.5, all inside the `gems` conda env.

## Global Constraints

- **Environment:** every command must run inside the `gems` conda env. The default `python3` on PATH is python 3.10 from a *different user's* anaconda (`/home/jyryu/anaconda3`) and has no access to this project's dependencies; `Rscript` is not on the default PATH at all. Always `conda activate gems` first, or use the absolute paths `/data2/home/farh89/miniconda3/envs/gems/bin/python` and `/data2/home/farh89/miniconda3/envs/gems/bin/Rscript`.
- **Python 3.6 / pandas 0.23.4.** f-strings are fine (the codebase uses them). `DataFrame.to_numpy()`, `df.iteritems()` replacements, and other pandas ≥0.24 APIs are NOT available — use `.values`.
- **This repo has no test framework, by deliberate decision** (commits `29afe18`, `a59e3bb` removed `Makefile`, `tox.ini`, and `tests/`). Do NOT create a `tests/` directory or add pytest. Verification in this plan uses inline heredoc scripts run against the four already-completed IBDMDB runs, which serve as real fixtures. The TDD cycle is preserved in spirit: every task writes its check first, runs it to see it fail, implements, then runs it to see it pass.
- **Fixture run used throughout:** `output_CD_ileum_Human2_tINIT_Human2_tINIT` — 36 patients, 2,285 genes in the union, 2,041 common to all 36. `PKM` = `ENSG00000067225` (present in 36/36), `SLC33A1` = `ENSG00000169359` (36/36).
- **Additive only.** When `--ko_targets` is absent, the pipeline must behave exactly as before. Do not modify `moma_project.R`, `moma_viz.R`, or any scoring code.
- **Display only.** This step computes no statistic and asserts no significance. Every figure carries the `[P0]` hypothesis-only caveat.
- **Style must match `moma_viz.R`:** Control `#009E73`, toward-control arrows `#D55E00` (vermilion), away arrows `grey70`, centre cross `#333333`, `theme_pub()` at `base_family = "Liberation Sans"`, saved at `width = 7.2, height = 6.2, dpi = 300, bg = "white"` as both `.png` and `.svg`.

---

### Task 1: Gene-name resolution helper

Resolves user-supplied gene names (HGNC symbols or model IDs) to model IDs using the run's own `gene_map.csv`, so the same command works against Recon2M.2 (Entrez) and Human2 (ENSG) runs.

**Files:**
- Modify: `tf_metabolism/utils.py` (add `_resolve_gene_targets` next to the existing `_effective_gene_map` helper)

**Interfaces:**
- Consumes: existing `_effective_gene_map(output_viz_dir) -> str` (already in `utils.py`)
- Produces: `_resolve_gene_targets(targets, gene_map_file) -> list of (model_id, symbol) tuples`
  - `targets`: list of str, each an HGNC symbol or a model gene ID
  - returns one tuple per resolved match, in input order; a symbol matching several model IDs yields several tuples
  - raises `ValueError` listing candidates when a name resolves to nothing

- [ ] **Step 1: Write the failing check**

Save as `/tmp/check_task1.py`:

```python
import sys
sys.path.insert(0, '/data2/home/farh89/projects/ibdmdb_gem/metabolic_modeling')
from tf_metabolism.utils import _resolve_gene_targets

GM = 'output_CD_ileum_Human2_tINIT_Human2_tINIT/viz/gene_map.csv'

# symbol -> model id
assert _resolve_gene_targets(['PKM'], GM) == [('ENSG00000067225', 'PKM')]
# model id passes through, symbol filled in
assert _resolve_gene_targets(['ENSG00000067225'], GM) == [('ENSG00000067225', 'PKM')]
# case-insensitive symbol
assert _resolve_gene_targets(['pkm'], GM) == [('ENSG00000067225', 'PKM')]
# order preserved across several inputs
out = _resolve_gene_targets(['PKM', 'SLC33A1'], GM)
assert [g for g, _ in out] == ['ENSG00000067225', 'ENSG00000169359'], out
# unresolvable -> ValueError naming the offender
try:
    _resolve_gene_targets(['NOTAGENE'], GM)
except ValueError as e:
    assert 'NOTAGENE' in str(e), str(e)
else:
    raise AssertionError('expected ValueError for NOTAGENE')

print('TASK 1 OK')
```

- [ ] **Step 2: Run it to verify it fails**

```bash
conda activate gems
cd /data2/home/farh89/projects/ibdmdb_gem/metabolic_modeling
python /tmp/check_task1.py
```

Expected: `ImportError: cannot import name '_resolve_gene_targets'`

- [ ] **Step 3: Implement**

Add to `tf_metabolism/utils.py`, immediately after `_effective_gene_map`:

```python
def _resolve_gene_targets(targets, gene_map_file):
    """Resolve user-supplied gene names to model gene IDs.

    Accepts either an HGNC symbol ('PKM') or a model gene ID
    ('ENSG00000067225' for Human-GEM runs, an Entrez ID for Recon2M.2 runs) and
    returns ``[(model_id, symbol), ...]`` in input order. A symbol mapping to
    several model IDs yields one tuple per match.

    NOTE: gene_map.csv's first column is named 'entrez' for historical reasons
    but holds ENSG IDs in Human-GEM runs — it is simply the model's gene ID.

    Raises ValueError naming the offending input (with near-miss candidates)
    when a name resolves to nothing, so a typo fails loud instead of silently
    producing no figure.
    """
    if not os.path.exists(gene_map_file):
        logging.warning('No gene map at %s — treating targets as model gene IDs.',
                        gene_map_file)
        return [(str(t), str(t)) for t in targets]

    gm = pd.read_csv(gene_map_file, dtype=str).fillna('')
    id_col, sym_col = gm.columns[0], gm.columns[1]
    ids = gm[id_col].tolist()
    syms = gm[sym_col].tolist()

    id_to_sym = dict(zip(ids, syms))
    sym_to_ids = {}
    for gid, sym in zip(ids, syms):
        if sym:
            sym_to_ids.setdefault(sym.upper(), []).append(gid)

    resolved = []
    for raw in targets:
        name = str(raw).strip()
        if not name:
            continue
        if name in id_to_sym:
            resolved.append((name, id_to_sym[name] or name))
            continue
        hits = sym_to_ids.get(name.upper())
        if hits:
            for gid in hits:
                resolved.append((gid, id_to_sym.get(gid) or name))
            continue
        near = [s for s in sym_to_ids if s.startswith(name.upper()[:3])][:10]
        raise ValueError(
            "Gene target '%s' not found in %s. Did you mean: %s"
            % (name, gene_map_file, ', '.join(sorted(near)) or '(no close matches)'))
    return resolved
```

- [ ] **Step 4: Run it to verify it passes**

```bash
python /tmp/check_task1.py
```

Expected: `TASK 1 OK`

- [ ] **Step 5: Commit**

```bash
git add tf_metabolism/utils.py
git commit -m "Add gene-name resolution helper for KO target figures"
```

---

### Task 2: The R renderer

Renders one figure per target gene from the per-patient projection CSVs. Written and verified standalone, before any Python calls it.

**Files:**
- Create: `tf_metabolism/r_based_viz/moma_cohort_viz.R`

**Interfaces:**
- Consumes: nothing from Task 1 (invoked directly by `Rscript` in this task)
- Produces: CLI contract used by Task 3 —
  `Rscript moma_cohort_viz.R --output-dir <run> --output-viz-dir <run>/viz --mode df --targets <comma-separated MODEL IDs> [--gene-map <path>]`
  Writes `<viz>/ko_cohort_<SYMBOL>.png` and `.svg`. Exits non-zero on fail-loud conditions.

- [ ] **Step 1: Write the failing check**

Save as `/tmp/check_task2.sh`:

```bash
#!/bin/bash
set -u
cd /data2/home/farh89/projects/ibdmdb_gem/metabolic_modeling
RUN=output_CD_ileum_Human2_tINIT_Human2_tINIT
R=tf_metabolism/r_based_viz/moma_cohort_viz.R
rm -f $RUN/viz/ko_cohort_PKM.* $RUN/viz/ko_cohort_SLC33A1.*

# 1. renders both requested genes
Rscript $R --output-dir $RUN --output-viz-dir $RUN/viz --mode df \
    --targets ENSG00000067225,ENSG00000169359 --gene-map $RUN/viz/gene_map.csv || exit 1
[ -f $RUN/viz/ko_cohort_PKM.png ] || { echo "FAIL: no PKM png"; exit 1; }
[ -f $RUN/viz/ko_cohort_PKM.svg ] || { echo "FAIL: no PKM svg"; exit 1; }
[ -f $RUN/viz/ko_cohort_SLC33A1.png ] || { echo "FAIL: no SLC33A1 png"; exit 1; }

# 2. a gene absent from every patient warns and skips, exit 0, no file
Rscript $R --output-dir $RUN --output-viz-dir $RUN/viz --mode df \
    --targets ENSGDOESNOTEXIST --gene-map $RUN/viz/gene_map.csv || \
    { echo "FAIL: absent gene should not be fatal"; exit 1; }
[ -f $RUN/viz/ko_cohort_ENSGDOESNOTEXIST.png ] && { echo "FAIL: wrote file for absent gene"; exit 1; }

# 3. missing reference file is fatal
Rscript $R --output-dir $RUN --output-viz-dir /tmp/definitely_empty_dir --mode df \
    --targets ENSG00000067225 2>/dev/null && { echo "FAIL: should have exited non-zero"; exit 1; }

echo "TASK 2 OK"
```

- [ ] **Step 2: Run it to verify it fails**

```bash
conda activate gems
mkdir -p /tmp/definitely_empty_dir
bash /tmp/check_task2.sh
```

Expected: fails at check 1 — `Rscript: ... moma_cohort_viz.R: No such file or directory`

- [ ] **Step 3: Implement**

Create `tf_metabolism/r_based_viz/moma_cohort_viz.R`:

```r
#!/usr/bin/env Rscript
# =============================================================================
# moma_cohort_viz.R — per-TARGET cohort view of MOMA knockouts
#
# Inverts moma_viz.R: instead of one figure per patient showing many genes,
# this draws one figure per GENE showing every patient's wild-type -> knockout
# shift against the nonIBD control cloud.
#
# All coordinates are re-sliced from the per-patient UMAP_Results_*.csv files
# written by moma_project.R — nothing is re-projected, so this figure shares
# the predict frame with the per-patient figures by construction.
#
# Usage:
#   Rscript moma_cohort_viz.R \
#       --output-dir     <run folder>     # for cohort_target_ranking.csv
#       --output-viz-dir <run folder>/viz # inputs and outputs
#       --mode           df               # default: df
#       --targets        ID1,ID2          # RESOLVED model gene IDs
#       [--gene-map      <viz>/gene_map.csv]   # labels only, no lookup of user input
#
# Outputs (per target):
#   ko_cohort_<SYMBOL>.png / .svg
#
# DISPLAY ONLY. Computes no statistic and asserts no significance; see
# FUTURE_REQUESTS [P0] (targets are hypothesis-only) and [P0b] (significance is
# never decided by 2D proximity).
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

set.seed(42)

option_list <- list(
  make_option("--output-dir",     type = "character",
              help = "run folder (for cohort_target_ranking.csv)"),
  make_option("--output-viz-dir", type = "character",
              help = "dir with UMAP_Results_*.csv and UMAP_predict_ref_<mode>.csv"),
  make_option("--mode",           type = "character", default = "df",
              help = "UMAP reference mode (default: df)"),
  make_option("--targets",        type = "character",
              help = "comma-separated RESOLVED model gene IDs"),
  make_option("--gene-map",       type = "character", default = NULL,
              help = "CSV with id,symbol columns — used only to label figures")
)

opt            <- parse_args(OptionParser(option_list = option_list))
output_dir     <- opt[["output-dir"]]
output_viz_dir <- opt[["output-viz-dir"]]
mode           <- opt[["mode"]]
targets        <- trimws(strsplit(opt[["targets"]], ",")[[1]])
gene_map_file  <- opt[["gene-map"]]

# ---------------------------------------------------------------------------
# Publication theme — identical to moma_viz.R so the two figure families match
# ---------------------------------------------------------------------------
theme_pub <- function(base_size = 10, base_family = "Liberation Sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.background   = element_rect(fill = "white"),
      panel.background  = element_rect(fill = "white"),
      axis.line         = element_line(color = "black", linewidth = 0.75),
      axis.ticks        = element_line(color = "black", linewidth = 0.75),
      axis.ticks.length = unit(6, "pt"),
      axis.title.x      = element_text(margin = margin(t = 5)),
      axis.title.y      = element_text(margin = margin(r = 5)),
      plot.title        = element_text(size = 12, face = "plain", hjust = 0.5),
      plot.subtitle     = element_text(size = 8,  hjust = 0.5, color = "grey25"),
      axis.title        = element_text(size = 12),
      axis.text         = element_text(size = 12, color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key        = element_rect(fill = "white", color = NA),
      legend.text       = element_text(size = 12),
      legend.title      = element_blank(),
      legend.position   = "right",
      plot.margin       = unit(c(0.15, 0.15, 0.15, 0.15), "cm")
    )
}

COLOR_MAP <- c(
  "Control"        = "#009E73",   # bluish green (Wong) — matches moma_viz.R
  "Patient WT"     = "#666666",
  "Patient KO"     = "#E69F00",
  "Control Center" = "#333333"
)

# ---------------------------------------------------------------------------
# Reference cloud (fail loud — the frame is undefined without it)
# ---------------------------------------------------------------------------
ref_file <- file.path(output_viz_dir, sprintf("UMAP_predict_ref_%s.csv", mode))
if (!file.exists(ref_file)) {
  stop("Predict-frame reference not found: ", ref_file,
       "\nRun the MOMA projection step first (moma_project.R).")
}
ref_df   <- read.csv(ref_file, row.names = 1, check.names = FALSE)
ctrl_df  <- ref_df[ref_df$condition == "Control", ]
center_x <- median(ctrl_df$UMAP1)
center_y <- median(ctrl_df$UMAP2)

# ---------------------------------------------------------------------------
# Per-patient projection results (fail loud if none)
# ---------------------------------------------------------------------------
res_files <- list.files(output_viz_dir,
                        pattern = "^UMAP_Results_MOMA_target_results_.*\\.csv$",
                        full.names = TRUE)
if (length(res_files) == 0) {
  stop("No UMAP_Results_MOMA_target_results_*.csv in ", output_viz_dir,
       "\nRun the MOMA projection step first (moma_project.R).")
}

per_patient <- list()
for (f in res_files) {
  sample_name <- sub("^UMAP_Results_MOMA_target_results_", "",
                     sub("\\.csv$", "", basename(f)))
  per_patient[[sample_name]] <- read.csv(f, row.names = 1, check.names = FALSE)
}

# ---------------------------------------------------------------------------
# Optional label + ranking-stat lookups
# ---------------------------------------------------------------------------
id2sym <- NULL
if (!is.null(gene_map_file) && file.exists(gene_map_file)) {
  gm <- read.csv(gene_map_file, colClasses = "character")
  id2sym <- setNames(gm[[2]], gm[[1]])
}

rank_df <- NULL
rank_file <- file.path(output_dir, "cohort_target_ranking.csv")
if (file.exists(rank_file)) {
  rank_df <- read.csv(rank_file, colClasses = c(Gene = "character"))
}

# ---------------------------------------------------------------------------
# Figure label per target. When two requested targets share an HGNC symbol (one
# symbol mapping to several model IDs), append the model ID so the output files
# stay distinct instead of silently overwriting each other.
# ---------------------------------------------------------------------------
sym_of <- function(gid) {
  if (!is.null(id2sym) && !is.na(id2sym[gid]) && nzchar(id2sym[gid]))
    id2sym[[gid]] else gid
}
label_of <- setNames(vapply(targets, sym_of, character(1)), targets)
for (d in names(which(table(label_of) > 1))) {
  idx <- which(label_of == d)
  label_of[idx] <- paste0(label_of[idx], "_", names(label_of)[idx])
}

# ---------------------------------------------------------------------------
# One figure per target
# ---------------------------------------------------------------------------
for (gene_id in targets) {

  rows <- list()
  n_missing <- 0
  for (sample_name in names(per_patient)) {
    d <- per_patient[[sample_name]]
    if (!(gene_id %in% rownames(d))) { n_missing <- n_missing + 1; next }
    r <- d[gene_id, ]
    if (is.na(r$wt_x) || is.na(r$X)) { n_missing <- n_missing + 1; next }
    rows[[sample_name]] <- data.frame(
      sample = sample_name,
      wt_x = r$wt_x, wt_y = r$wt_y,
      ko_x = r$X,    ko_y = r$Y,
      # Toward/away is judged on the UMAP recovery axis — the SAME frame the
      # figure draws. Using recovery_proj_pca here would colour arrows by a
      # quantity the viewer cannot see, so arrows would visibly point toward
      # control while labelled "away". recovery_proj_pca remains the scoring
      # quantity everywhere else; this choice is local to the display.
      toward = !is.na(r$recovery_proj_umap) && r$recovery_proj_umap > 0,
      stringsAsFactors = FALSE)
  }

  symbol <- label_of[[gene_id]]

  if (length(rows) == 0) {
    warning(sprintf("Target %s (%s) is absent from all %d patient models — skipping.",
                    symbol, gene_id, length(per_patient)), immediate. = TRUE)
    next
  }

  arrow_df <- do.call(rbind, rows)
  n_shown  <- nrow(arrow_df)
  n_toward <- sum(arrow_df$toward)

  stat_txt <- ""
  if (!is.null(rank_df) && gene_id %in% rank_df$Gene) {
    rr <- rank_df[rank_df$Gene == gene_id, ][1, ]
    stat_txt <- sprintf("  |  rra_fdr = %.3g, recurrence = %.2f",
                        rr$rra_fdr, rr$recurrence_prop)
  }
  miss_txt <- if (n_missing > 0)
    sprintf("  |  %d patient(s) dropped (gene not in their model)", n_missing) else ""

  subtitle <- sprintf(
    "%d/%d patients moved toward control (UMAP recovery axis)%s%s\nDisplay only — MOMA targets are hypothesis-only (FUTURE_REQUESTS [P0]); significance is never read off a 2D embedding.",
    n_toward, n_shown, stat_txt, miss_txt)

  cross_df <- data.frame(x = center_x, y = center_y, group = "Control Center")

  p <- ggplot() +
    stat_ellipse(data = ctrl_df,
                 aes(x = UMAP1, y = UMAP2, color = "Control", fill = "Control"),
                 type = "norm", level = 0.80, linetype = "dashed", linewidth = 0.3,
                 geom = "polygon", alpha = 0.06, show.legend = FALSE) +
    geom_point(data = ctrl_df,
               aes(x = UMAP1, y = UMAP2, color = "Control"),
               alpha = 0.7, size = 2.4, shape = 16) +
    geom_segment(data = arrow_df[!arrow_df$toward, , drop = FALSE],
                 aes(x = wt_x, y = wt_y, xend = ko_x, yend = ko_y),
                 color = "grey70", alpha = 0.55, linewidth = 0.35,
                 arrow = arrow(length = unit(0.12, "cm"), type = "open"),
                 show.legend = FALSE) +
    geom_segment(data = arrow_df[arrow_df$toward, , drop = FALSE],
                 aes(x = wt_x, y = wt_y, xend = ko_x, yend = ko_y),
                 color = "#D55E00", alpha = 0.65, linewidth = 0.35,
                 arrow = arrow(length = unit(0.12, "cm"), type = "open"),
                 show.legend = FALSE) +
    geom_point(data = arrow_df, aes(x = wt_x, y = wt_y, color = "Patient WT"),
               size = 2.0, shape = 1, stroke = 0.9) +
    geom_point(data = arrow_df, aes(x = ko_x, y = ko_y, color = "Patient KO"),
               size = 2.2, shape = 16, alpha = 0.9) +
    geom_point(data = cross_df, aes(x = x, y = y, color = group),
               shape = 3, size = 3.5, stroke = 0.9) +
    scale_color_manual(values = COLOR_MAP, name = NULL) +
    scale_fill_manual(values = COLOR_MAP, guide = "none") +
    labs(x = "UMAP 1", y = "UMAP 2",
         title = bquote("Cohort response to " * italic(.(symbol)) * " knockout"),
         subtitle = subtitle) +
    theme_pub() +
    theme(aspect.ratio = 1)

  safe <- gsub("[^A-Za-z0-9_.-]", "_", symbol)
  out_png <- file.path(output_viz_dir, sprintf("ko_cohort_%s.png", safe))
  out_svg <- file.path(output_viz_dir, sprintf("ko_cohort_%s.svg", safe))
  ggsave(out_png, p, width = 7.2, height = 6.2, dpi = 300, bg = "white")
  ggsave(out_svg, p, width = 7.2, height = 6.2, bg = "white")
  message(sprintf("[MOMA-Cohort] %s: %d/%d toward control, %d dropped -> %s",
                  symbol, n_toward, n_shown, n_missing, basename(out_png)))
}
```

- [ ] **Step 4: Run it to verify it passes**

```bash
bash /tmp/check_task2.sh
```

Expected: `TASK 2 OK`

- [ ] **Step 5: Verify the frame matches the per-patient figures**

This is the correctness check that matters — the new figure must sit in the same predict frame. Run:

```bash
python - <<'EOF'
import glob, os
import pandas as pd
viz = 'output_CD_ileum_Human2_tINIT_Human2_tINIT/viz'
g = 'ENSG00000067225'   # PKM
n = 0
for f in sorted(glob.glob(os.path.join(viz, 'UMAP_Results_MOMA_target_results_*.csv'))):
    d = pd.read_csv(f, index_col=0)
    assert g in d.index, os.path.basename(f)
    r = d.loc[g]
    assert pd.notnull(r['wt_x']) and pd.notnull(r['X'])
    n += 1
print('PKM present with WT+KO coords in %d/%d patients' % (n, 36))
assert n == 36
print('FRAME CHECK OK')
EOF
```

Expected: `PKM present with WT+KO coords in 36/36 patients` then `FRAME CHECK OK`. Open `ko_cohort_PKM.png` and confirm 36 arrows, an open WT dot and filled KO dot at each end, the green control cloud, and the `[P0]` caveat in the subtitle.

- [ ] **Step 6: Commit**

```bash
git add tf_metabolism/r_based_viz/moma_cohort_viz.R
git commit -m "Add per-target MOMA cohort view renderer"
```

---

### Task 3: Python wrapper

Resolves names via Task 1 and shells out to Task 2's script, mirroring the R branch of `visualize_moma_umap`.

**Files:**
- Modify: `tf_metabolism/utils.py` (add `visualize_cohort_ko_umap` after `visualize_moma_umap`)

**Interfaces:**
- Consumes: `_resolve_gene_targets(targets, gene_map_file)` (Task 1), `_effective_gene_map(output_viz_dir)` and `_find_rscript()` and `_r_script(name)` (all pre-existing in `utils.py`), and the CLI contract of `moma_cohort_viz.R` (Task 2)
- Produces: `visualize_cohort_ko_umap(output_dir, output_viz_dir, targets, mode='df') -> None`, used by Task 4

- [ ] **Step 1: Write the failing check**

Save as `/tmp/check_task3.py`:

```python
import os, sys, glob
sys.path.insert(0, '/data2/home/farh89/projects/ibdmdb_gem/metabolic_modeling')
from tf_metabolism.utils import visualize_cohort_ko_umap

RUN = 'output_CD_ileum_Human2_tINIT_Human2_tINIT'
VIZ = os.path.join(RUN, 'viz')
for f in glob.glob(os.path.join(VIZ, 'ko_cohort_*')):
    os.remove(f)

# symbols in, figures out
visualize_cohort_ko_umap(RUN, VIZ, ['PKM', 'SLC33A1'], mode='df')
assert os.path.exists(os.path.join(VIZ, 'ko_cohort_PKM.png')), 'no PKM png'
assert os.path.exists(os.path.join(VIZ, 'ko_cohort_PKM.svg')), 'no PKM svg'
assert os.path.exists(os.path.join(VIZ, 'ko_cohort_SLC33A1.png')), 'no SLC33A1 png'

# a bad name fails loud before R is ever invoked
try:
    visualize_cohort_ko_umap(RUN, VIZ, ['NOTAGENE'], mode='df')
except ValueError as e:
    assert 'NOTAGENE' in str(e), str(e)
else:
    raise AssertionError('expected ValueError for NOTAGENE')

# empty target list is a no-op, not an error
visualize_cohort_ko_umap(RUN, VIZ, [], mode='df')

print('TASK 3 OK')
```

- [ ] **Step 2: Run it to verify it fails**

```bash
conda activate gems
cd /data2/home/farh89/projects/ibdmdb_gem/metabolic_modeling
python /tmp/check_task3.py
```

Expected: `ImportError: cannot import name 'visualize_cohort_ko_umap'`

- [ ] **Step 3: Implement**

Add to `tf_metabolism/utils.py`, immediately after `visualize_moma_umap`:

```python
def visualize_cohort_ko_umap(output_dir, output_viz_dir, targets, mode='df'):
    """Render one cohort-level figure per knocked-out target gene.

    Inverts the per-patient MOMA view: for each gene in ``targets``, draws every
    condition2 patient's wild-type -> knockout shift against the nonIBD control
    cloud, re-slicing the per-patient ``UMAP_Results_*.csv`` files written by
    moma_project.R. Nothing is re-projected and no MOMA simulation is re-run, so
    the figure shares the predict frame with the per-patient figures.

    Parameters
    ----------
    output_dir     : run folder (read for cohort_target_ranking.csv, if present)
    output_viz_dir : viz folder holding the projection results; figures land here
    targets        : list of HGNC symbols and/or model gene IDs
    mode           : UMAP reference mode (default 'df')

    Writes ``{output_viz_dir}/ko_cohort_<SYMBOL>.png`` (+ .svg) per target.
    Raises ValueError on an unresolvable gene name and RuntimeError if the R
    renderer exits non-zero. Display only — see FUTURE_REQUESTS [P0].
    """
    if not targets:
        return

    gene_map = _effective_gene_map(output_viz_dir)
    resolved = _resolve_gene_targets(targets, gene_map)
    if not resolved:
        logging.warning('No usable KO targets after resolution — skipping.')
        return

    gene_ids = [gid for gid, _ in resolved]
    logging.info('[MOMA-Cohort] Rendering %d target(s): %s',
                 len(resolved), ', '.join('%s (%s)' % (s, g) for g, s in resolved))

    cmd = [
        _find_rscript(), _r_script('moma_cohort_viz.R'),
        '--output-dir',     output_dir,
        '--output-viz-dir', output_viz_dir,
        '--mode',           mode,
        '--targets',        ','.join(gene_ids),
    ]
    if os.path.exists(gene_map):
        cmd += ['--gene-map', gene_map]

    logging.info('[R] Running: %s', ' '.join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError('moma_cohort_viz.R exited with code %d' % result.returncode)
    logging.info('[R] MOMA cohort target visualization done.')
```

If `subprocess` is not already imported at module level in `utils.py`, add `import subprocess` to the imports at the top rather than inside the function (`visualize_moma_umap` imports it locally; module-level is cleaner and harmless).

- [ ] **Step 4: Run it to verify it passes**

```bash
python /tmp/check_task3.py
```

Expected: `TASK 3 OK`

- [ ] **Step 5: Commit**

```bash
git add tf_metabolism/utils.py
git commit -m "Add Python wrapper for per-target MOMA cohort view"
```

---

### Task 4: CLI flag and pipeline wiring

Exposes the step as `--ko_targets` and calls it from `run_pipeline`, guarded so an absent flag changes nothing.

**Files:**
- Modify: `tf_metabolism/utils.py` — `argument_parser`, after the `--target_nboot` line
- Modify: `tf_metabolism/__main__.py:241-245` (args tuple), `:354-356` (signature), `:570-574` (call site)

**Interfaces:**
- Consumes: `visualize_cohort_ko_umap(output_dir, output_viz_dir, targets, mode='df')` (Task 3)
- Produces: `--ko_targets` CLI flag; `run_pipeline(..., target_nboot, ko_targets, start)`

- [ ] **Step 1: Write the failing check**

Save as `/tmp/check_task4.py`:

```python
import sys, inspect
sys.path.insert(0, '/data2/home/farh89/projects/ibdmdb_gem/metabolic_modeling')
from tf_metabolism.utils import argument_parser
from tf_metabolism.__main__ import run_pipeline

p = argument_parser()

# flag exists, defaults to None (absent -> step skipped)
ns = p.parse_args(['-o', 'x', '-c1', 'a.csv', '-c2', 'b.csv'])
assert ns.ko_targets is None, ns.ko_targets

# space-separated
ns = p.parse_args(['-o', 'x', '-c1', 'a.csv', '-c2', 'b.csv',
                   '--ko_targets', 'PKM', 'ENO1'])
assert ns.ko_targets == ['PKM', 'ENO1'], ns.ko_targets

# comma-separated, same as --model/--mem_method convention
ns = p.parse_args(['-o', 'x', '-c1', 'a.csv', '-c2', 'b.csv',
                   '--ko_targets', 'PKM,ENO1'])
assert ns.ko_targets == ['PKM,ENO1'], ns.ko_targets   # split happens downstream

# run_pipeline accepts the new parameter, positioned before `start`
params = list(inspect.signature(run_pipeline).parameters)
assert 'ko_targets' in params, params
assert params.index('ko_targets') == params.index('start') - 1, params

print('TASK 4 OK')
```

- [ ] **Step 2: Run it to verify it fails**

```bash
conda activate gems
cd /data2/home/farh89/projects/ibdmdb_gem/metabolic_modeling
python /tmp/check_task4.py
```

Expected: `AttributeError: 'Namespace' object has no attribute 'ko_targets'`

- [ ] **Step 3a: Add the flag**

In `tf_metabolism/utils.py`, in `argument_parser`, immediately after the `--target_nboot` argument:

```python
    parser.add_argument('--ko_targets', nargs='+', default=None,
                        help="Gene(s) to render a cohort-level knockout figure for, "
                             "one or more (space- or comma-separated). Accepts HGNC "
                             "symbols or model gene IDs. Each produces "
                             "viz/ko_cohort_<SYMBOL>.png showing every patient's "
                             "WT->KO shift vs the control cloud. Omit to skip. "
                             "Display only: MOMA targets are hypothesis-only.")
```

- [ ] **Step 3b: Thread it through the call site**

In `tf_metabolism/__main__.py`, in the args tuple (currently ending `options.target_nperm, options.target_nboot, start)`):

```python
        args = (model_name, mem, combo_out, omics_file1, omics_file2,
                options.flux_pval, options.flux_filter, options.flux_rbc,
                options.target_fdr, options.target_effect_min, options.target_topk,
                options.target_method, options.target_switch_n,
                options.target_nperm, options.target_nboot,
                options.ko_targets, start)
```

- [ ] **Step 3c: Add it to the signature**

In `tf_metabolism/__main__.py`, change the `run_pipeline` signature line from
`target_nboot, start):` to:

```python
                 target_nboot, ko_targets, start):
```

- [ ] **Step 3d: Call it**

In `tf_metabolism/__main__.py`, inside the existing `if target_reactions_df:` block, after `visualize_moma_umap(...)`:

```python
    if target_reactions_df:
        run_moma_results_umap(output_dir, output_viz_dir, target_reactions_df, mode='df', use_r=True)
        visualize_moma_umap(output_viz_dir, mode='df', use_r=True)
        ## Per-target cohort KO view (opt-in via --ko_targets) → viz/ko_cohort_<GENE>.png
        ## Purely additive: absent flag = no-op. Display only, see FUTURE_REQUESTS [P0].
        if ko_targets:
            try:
                visualize_cohort_ko_umap(output_dir, output_viz_dir, ko_targets, mode='df')
            except Exception as e:
                logging.warning('Per-target cohort KO view skipped: %s', e)
    else:
        logging.warning('No DF reactions — skipping MOMA projection.')
```

Also add `visualize_cohort_ko_umap` to the `from tf_metabolism.utils import (...)` block at `__main__.py:24-29`.

Note the flag arrives as a list that may contain comma-joined entries (`['PKM,ENO1']`). Normalize inside `visualize_cohort_ko_umap` by replacing its first line with:

```python
    if not targets:
        return
    flat = []
    for t in targets:
        flat.extend([x for x in str(t).split(',') if x.strip()])
    targets = flat
```

- [ ] **Step 4: Run it to verify it passes**

```bash
python /tmp/check_task4.py
```

Expected: `TASK 4 OK`

- [ ] **Step 5: Verify comma-splitting end to end**

```bash
python - <<'EOF'
import os, sys, glob
sys.path.insert(0, '/data2/home/farh89/projects/ibdmdb_gem/metabolic_modeling')
from tf_metabolism.utils import visualize_cohort_ko_umap
RUN='output_CD_ileum_Human2_tINIT_Human2_tINIT'; VIZ=os.path.join(RUN,'viz')
for f in glob.glob(os.path.join(VIZ,'ko_cohort_*')): os.remove(f)
visualize_cohort_ko_umap(RUN, VIZ, ['PKM,SLC33A1'], mode='df')   # one comma-joined arg
assert os.path.exists(os.path.join(VIZ,'ko_cohort_PKM.png'))
assert os.path.exists(os.path.join(VIZ,'ko_cohort_SLC33A1.png'))
print('COMMA SPLIT OK')
EOF
```

Expected: `COMMA SPLIT OK`

- [ ] **Step 6: Verify the absent flag is a true no-op**

```bash
conda activate gems
cd /data2/home/farh89/projects/ibdmdb_gem/metabolic_modeling
python ./tf_metabolism.py --help | grep -A 6 "ko_targets"
git stash list > /dev/null   # no-op, keeps the shell honest
python -c "
import sys; sys.path.insert(0,'.')
from tf_metabolism.utils import argument_parser
ns = argument_parser().parse_args(['-o','x','-c1','a.csv','-c2','b.csv'])
print('ko_targets default =', ns.ko_targets)
assert ns.ko_targets is None
print('NO-OP DEFAULT OK')"
```

Expected: the help text for `--ko_targets`, then `NO-OP DEFAULT OK`.

- [ ] **Step 7: Commit**

```bash
git add tf_metabolism/utils.py tf_metabolism/__main__.py
git commit -m "Wire --ko_targets flag into pipeline for per-target cohort KO view"
```

---

### Task 5: Documentation

The repo's convention is that CLAUDE.md carries a numbered section per local extension, and WORKLOG.md gets a dated session entry.

**Files:**
- Modify: `CLAUDE.md` (new §3i, and a §8 pointer)
- Modify: `WORKLOG.md` (new dated entry at the top — newest first)
- Modify: `README.md` (document the new flag alongside the other CLI options)

**Interfaces:**
- Consumes: the finished feature from Tasks 1-4
- Produces: nothing consumed by later tasks

- [ ] **Step 1: Add CLAUDE.md §3i**

Insert after §3h, before the "### NOT yet implemented" heading:

```markdown
### 3i. Per-target MOMA cohort view (implemented, CORE pipeline step, opt-in)
`--ko_targets GENE[,GENE...]` renders `viz/ko_cohort_<SYMBOL>.png` (+ .svg) per gene:
every condition2 patient's WT→KO shift as an arrow over the nonIBD control cloud —
the INVERSE of the per-patient `UMAP_MOMA_target_results_<id>.png`. Accepts HGNC
symbols or model gene IDs (resolved in Python via `viz/gene_map.csv`, so the same
command works for Recon2M.2/Entrez and Human2/ENSG runs); an unresolvable name fails
loud. Renderer: `r_based_viz/moma_cohort_viz.R` (R/ggplot2 only, no matplotlib
fallback — deliberately avoids the dual-backend divergence of FUTURE_REQUESTS [P2]).
Wrapper: `utils.visualize_cohort_ko_umap`. Called from `run_pipeline` inside the
existing `if target_reactions_df:` block; **absent flag = complete no-op**.
Re-slices the per-patient `UMAP_Results_*.csv` by gene — nothing is re-projected,
`moma_project.R` is untouched, and NO MOMA/Gurobi re-run is needed, so it renders
for already-completed runs. Arrow colour (toward/away) uses `recovery_proj_umap`,
NOT `recovery_proj_pca`: the figure draws the UMAP frame, so colouring by an unseen
quantity would contradict the visible arrows. `recovery_proj_pca` remains the
scoring quantity everywhere else. Gene coverage varies per patient (tINIT builds a
context model each), so patients lacking the gene are dropped and counted in the
subtitle. DISPLAY ONLY — carries the [P0] hypothesis-only caveat.
```

- [ ] **Step 2: Add the CLAUDE.md §8 pointer**

Append to the §8 pointer list:

```markdown
- `tf_metabolism/r_based_viz/moma_cohort_viz.R` — per-target cohort KO figure (§3i);
  core pipeline step, opt-in via `--ko_targets`, renders from existing run outputs.
```

- [ ] **Step 3: Add the WORKLOG entry**

Insert directly below the `Session-by-session log...` header line, above the newest existing session:

```markdown
## Session: 2026-07-29b (per-target MOMA cohort view — new core pipeline step)

Added `--ko_targets`: one figure per knocked-out gene showing every condition2
patient's WT→KO shift against the control cloud, inverting the per-patient MOMA view.
Design spec: `docs/superpowers/specs/2026-07-29-per-target-moma-cohort-view-design.md`;
plan: `docs/superpowers/plans/2026-07-29-per-target-moma-cohort-view.md`.

- **No re-simulation.** `moma_project.R` already projects EVERY gene for EVERY patient;
  the new step re-slices `viz/UMAP_Results_*.csv` by gene instead of by patient. So it
  renders for the four completed IBDMDB runs with no Gurobi re-run, and `moma_project.R`
  (and the [P0b] fix in it) is untouched.
- **Opt-in.** Flag absent = complete no-op; existing runs behave identically.
- **R/ggplot2 only** — deliberately no matplotlib fallback, to avoid recreating the
  dual-backend silent divergence documented as [P2].
- **Arrow toward/away uses `recovery_proj_umap`, not `recovery_proj_pca`** — the figure
  draws the UMAP frame, so colouring by the PCA score would produce arrows visibly
  pointing toward control while labelled "away". Display-local; scoring unchanged.
- **Gene coverage varies per patient** (tINIT context models): CD-ileum Human2 has 2,041
  genes common to all 36 patients but 2,285 in the union. Patients lacking the target are
  dropped and the count is reported in the subtitle.
- Display only; subtitle carries the [P0] hypothesis-only caveat.
```

- [ ] **Step 4: Document the flag in README.md**

Add to the CLI options section, matching the surrounding format:

```markdown
`--ko_targets GENE [GENE ...]` — render a cohort-level knockout figure per gene
(`viz/ko_cohort_<SYMBOL>.png`): every condition2 patient's wild-type → knockout shift
against the nonIBD control cloud. Accepts HGNC symbols or model gene IDs, space- or
comma-separated. Omit to skip. Requires a completed MOMA projection step; renders from
existing outputs without re-running MOMA. Display only — MOMA targets are
hypothesis-only (see FUTURE_REQUESTS [P0]).
```

- [ ] **Step 5: Verify the docs match the code**

```bash
cd /data2/home/farh89/projects/ibdmdb_gem/metabolic_modeling
grep -n "ko_targets" CLAUDE.md README.md WORKLOG.md tf_metabolism/utils.py tf_metabolism/__main__.py
grep -n "3i" CLAUDE.md | head -3
```

Expected: `--ko_targets` appears in all three docs and both source files, and §3i exists in CLAUDE.md.

- [ ] **Step 6: Commit**

```bash
git add CLAUDE.md WORKLOG.md README.md
git commit -m "Document per-target MOMA cohort view (--ko_targets)"
```

---

## Final verification

Run after all tasks, from a clean checkout of the branch:

```bash
conda activate gems
cd /data2/home/farh89/projects/ibdmdb_gem/metabolic_modeling
RUN=output_CD_ileum_Human2_tINIT_Human2_tINIT

# 1. renders for both cohorts, including the one shared CD/UC target
python -c "
import sys; sys.path.insert(0,'.')
from tf_metabolism.utils import visualize_cohort_ko_umap
for run in ['output_CD_ileum_Human2_tINIT_Human2_tINIT',
            'output_UC_rectum_Human2_tINIT_Human2_tINIT']:
    visualize_cohort_ko_umap(run, run+'/viz', ['SLC33A1'], mode='df')
print('BOTH COHORTS OK')"

# 2. arrow count equals patient count (36 CD / 20 UC), verified from source data
python - <<'EOF'
import glob, os
import pandas as pd
for run, n_expected in [('output_CD_ileum_Human2_tINIT_Human2_tINIT', 36),
                        ('output_UC_rectum_Human2_tINIT_Human2_tINIT', 20)]:
    fs = sorted(glob.glob(os.path.join(run, 'viz', 'UMAP_Results_MOMA_target_results_*.csv')))
    n = sum('ENSG00000169359' in pd.read_csv(f, index_col=0).index for f in fs)
    print('%s: SLC33A1 in %d/%d patients' % (run, n, len(fs)))
    assert len(fs) == n_expected, (len(fs), n_expected)
print('ARROW COUNT SOURCE OK')
EOF

# 3. the four PNGs exist and are non-empty
ls -l output_*_Human2_tINIT_Human2_tINIT/viz/ko_cohort_*.png
```

Then open `ko_cohort_SLC33A1.png` for both cohorts and confirm visually: green control cloud with a dark `+` at its centroid, one arrow per patient, vermilion arrows toward control and grey away, and the subtitle showing the toward/total count plus the `[P0]` caveat.

## Notes for the implementer

- **Do not** create a `tests/` directory or add pytest — their removal was deliberate (`29afe18`, `a59e3bb`).
- **Do not** touch `moma_project.R`, `moma_viz.R`, or anything in `cohort_target_ranking.py`. If a change there looks necessary, stop and raise it — it means the re-slice approach has a gap the spec missed.
- The `output_*/` folders are **gitignored and not recoverable from git**. The checks above delete only `ko_cohort_*` files, which are regenerable; don't widen those globs.
- `gene_map.csv`'s first column is named `entrez` but holds ENSG IDs in Human-GEM runs. This is a known upstream misnomer — read it positionally, as Task 1 does.
