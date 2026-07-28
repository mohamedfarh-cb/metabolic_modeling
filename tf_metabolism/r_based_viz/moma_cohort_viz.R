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

# ---------------------------------------------------------------------------
# Legend / colour key.
#
# Arrow COLOUR is the direction encoding (toward vs away from the control
# cloud), so it must be a LABELLED legend entry — not a bare `color=` argument
# outside aes(). Both arrow levels therefore go through the same colour scale
# as the points, and `override.aes` below gives each key the right glyph:
# points draw as points (linetype 0), arrows draw as lines (shape NA). Without
# that split every key would show a dot AND a line and the encoding would stay
# unreadable.
#
# "Patient KO" was #E69F00 (orange), which sat next to the #D55E00 (vermillion)
# toward-control arrow and invited the reader to map orange-line -> KO. It is
# now #0072B2 (blue, Wong) so no point colour is adjacent to an arrow colour.
# Arrow colours themselves are unchanged. "Control" #009E73 is unchanged.
# ---------------------------------------------------------------------------
LVL_CTRL   <- "Control"
LVL_CENTER <- "Control Center"
LVL_WT     <- "Patient WT"
LVL_KO     <- "Patient KO"
LVL_TOWARD <- "Arrow: KO toward control"
LVL_AWAY   <- "Arrow: KO away from control"

COLOR_MAP <- c(
  "#009E73",   # Control        — bluish green (Wong), matches moma_viz.R
  "#333333",   # Control Center
  "#666666",   # Patient WT     — open dot
  "#0072B2",   # Patient KO     — filled dot (was #E69F00; collided with arrows)
  "#D55E00",   # toward-control arrow (unchanged)
  "grey70"     # away-from-control arrow (unchanged)
)
names(COLOR_MAP) <- c(LVL_CTRL, LVL_CENTER, LVL_WT, LVL_KO, LVL_TOWARD, LVL_AWAY)

# Per-level key glyph: shape NA + linetype 1 => a LINE key (arrows);
#                      shape n  + linetype 0 => a POINT key (dots/cross).
KEY_SHAPE    <- c(16, 3, 1, 16, NA, NA)
KEY_LINETYPE <- c(0,  0, 0, 0,  1,  1)
names(KEY_SHAPE) <- names(KEY_LINETYPE) <- names(COLOR_MAP)

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

# A legacy/partial cohort_target_ranking.csv (no rra_fdr / recurrence_prop, or no
# Gene column at all) must NOT abort the run: sprintf() on a NULL column errors and
# would kill every remaining target in the same invocation. Validate up front and
# fall back to "no stats in the subtitle", which the spec already allows for a
# missing file.
RANK_COLS <- c("Gene", "rra_fdr", "recurrence_prop")
rank_df   <- NULL
rank_file <- file.path(output_dir, "cohort_target_ranking.csv")
if (file.exists(rank_file)) {
  rank_df <- tryCatch(read.csv(rank_file, colClasses = c(Gene = "character")),
                      error = function(e) {
                        warning(sprintf("Could not read %s (%s) — subtitle stats omitted.",
                                        rank_file, conditionMessage(e)), immediate. = TRUE)
                        NULL
                      })
  if (!is.null(rank_df) && !all(RANK_COLS %in% names(rank_df))) {
    warning(sprintf("%s lacks column(s) %s — subtitle stats omitted.", rank_file,
                    paste(setdiff(RANK_COLS, names(rank_df)), collapse = ", ")),
            immediate. = TRUE)
    rank_df <- NULL
  }
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
    stat_txt <- sprintf(" | rra_fdr = %.3g, recurrence = %.2f",
                        rr$rra_fdr, rr$recurrence_prop)
  }
  miss_txt <- if (n_missing > 0)
    sprintf(" | %d patient(s) dropped (gene not in their model)", n_missing) else ""

  # The stats line has no fixed width: symbol length + %.3g values as long as
  # "1.08e-110" + a two-digit dropped count. With plot.title.position = "plot"
  # the subtitle is centred on the FULL device width, so an over-long line clips
  # at BOTH ends (observed: ST6GALNAC3 lost its leading "28/35" and trailing ")").
  # Wrap it. SUBTITLE_WRAP is set below the fixed [P0] caveat line's own length
  # (123 chars, which renders intact at 7.2in / size 8), so both always fit.
  SUBTITLE_WRAP <- 100
  line1 <- sprintf("%d/%d patients moved toward control (UMAP recovery axis)%s%s",
                   n_toward, n_shown, stat_txt, miss_txt)
  line1 <- paste(strwrap(line1, width = SUBTITLE_WRAP), collapse = "\n")

  subtitle <- paste0(
    line1,
    "\nDisplay only — MOMA targets are hypothesis-only (FUTURE_REQUESTS [P0]); significance is never read off a 2D embedding.")

  cross_df <- data.frame(x = center_x, y = center_y, group = LVL_CENTER)

  # Arrow direction is mapped THROUGH aes() so it becomes a named legend entry.
  arrow_df$direction <- ifelse(arrow_df$toward, LVL_TOWARD, LVL_AWAY)

  # Only the levels actually drawn get a key (a cohort can be 100% toward), so
  # the override.aes vectors are subset by name to stay length-matched.
  shown_levels <- names(COLOR_MAP)[names(COLOR_MAP) %in%
                                   c(LVL_CTRL, LVL_CENTER, LVL_WT, LVL_KO,
                                     unique(arrow_df$direction))]

  p <- ggplot() +
    stat_ellipse(data = ctrl_df,
                 aes(x = UMAP1, y = UMAP2, color = LVL_CTRL, fill = LVL_CTRL),
                 type = "norm", level = 0.80, linetype = "dashed", linewidth = 0.3,
                 geom = "polygon", alpha = 0.06, show.legend = FALSE) +
    geom_point(data = ctrl_df,
               aes(x = UMAP1, y = UMAP2, color = LVL_CTRL),
               alpha = 0.7, size = 2.4, shape = 16) +
    geom_segment(data = arrow_df[!arrow_df$toward, , drop = FALSE],
                 aes(x = wt_x, y = wt_y, xend = ko_x, yend = ko_y, color = direction),
                 alpha = 0.55, linewidth = 0.35,
                 arrow = arrow(length = unit(0.12, "cm"), type = "open")) +
    geom_segment(data = arrow_df[arrow_df$toward, , drop = FALSE],
                 aes(x = wt_x, y = wt_y, xend = ko_x, yend = ko_y, color = direction),
                 alpha = 0.65, linewidth = 0.35,
                 arrow = arrow(length = unit(0.12, "cm"), type = "open")) +
    geom_point(data = arrow_df, aes(x = wt_x, y = wt_y, color = LVL_WT),
               size = 2.0, shape = 1, stroke = 0.9) +
    geom_point(data = arrow_df, aes(x = ko_x, y = ko_y, color = LVL_KO),
               size = 2.2, shape = 16, alpha = 0.9) +
    geom_point(data = cross_df, aes(x = x, y = y, color = group),
               shape = 3, size = 3.5, stroke = 0.9) +
    scale_color_manual(
      values = COLOR_MAP, name = NULL, breaks = shown_levels,
      guide  = guide_legend(
        override.aes = list(shape    = unname(KEY_SHAPE[shown_levels]),
                            linetype = unname(KEY_LINETYPE[shown_levels]),
                            size     = 2.6, linewidth = 0.9, alpha = 1))) +
    scale_fill_manual(values = COLOR_MAP, guide = "none") +
    labs(x = "UMAP 1", y = "UMAP 2",
         title = bquote("Cohort response to " * italic(.(symbol)) * " knockout"),
         subtitle = subtitle) +
    theme_pub() +
    theme(aspect.ratio = 1,
          # Anchor title/subtitle to the FULL plot width (device), not just the
          # panel: with a legend on the right, "panel"-anchored centering (the
          # ggplot2 default) leaves too little room for the wide [P0] caveat
          # line and clips its left edge off-canvas. "plot" anchoring gives the
          # subtitle the full device width to center in, so the fixed caveat
          # sentence renders fully on-canvas without being shortened.
          plot.title.position = "plot")

  safe <- gsub("[^A-Za-z0-9_.-]", "_", symbol)
  out_png <- file.path(output_viz_dir, sprintf("ko_cohort_%s.png", safe))
  out_svg <- file.path(output_viz_dir, sprintf("ko_cohort_%s.svg", safe))
  ggsave(out_png, p, width = 7.2, height = 6.2, dpi = 300, bg = "white")
  ggsave(out_svg, p, width = 7.2, height = 6.2, bg = "white")
  message(sprintf("[MOMA-Cohort] %s: %d/%d toward control, %d dropped -> %s",
                  symbol, n_toward, n_shown, n_missing, basename(out_png)))
}
