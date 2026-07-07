#!/usr/bin/env Rscript
# =============================================================================
# moma_viz.R — MOMA simulation result visualization
#
# Reads UMAP_Results_*.csv produced by run_moma_results_umap() (Python)
# and generates publication-quality plots with ggplot2 + ggrepel.
#
# Usage:
#   Rscript moma_viz.R \
#       --output-viz-dir <path>    # dir containing UMAP_Results_*.csv
#       [--mode          dc]       # UMAP reference mode (default: dc)
#
# Outputs (per sample):
#   UMAP_{basename}.png / .svg    — scatter plot
# Shared output:
#   target_genes.csv              — top-10 target genes per sample
#                                   (ranked by recovery_proj_pca toward Control;
#                                    columns: Gene, Sample, rank, recovery_proj_pca,
#                                    delta_distance, distance_to_center)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(ggrepel)
})

set.seed(42)

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--output-viz-dir", type = "character",
              help = "directory with UMAP_Results_*.csv and UMAP_{mode}_raw_data.csv"),
  make_option("--mode", type = "character", default = "dc",
              help = "UMAP reference mode (default: dc)"),
  make_option("--gene-map", type = "character", default = NULL,
              help = "CSV with 'entrez,symbol' columns to label genes by HGNC symbol")
)

opt            <- parse_args(OptionParser(option_list = option_list))
output_viz_dir <- opt[["output-viz-dir"]]
mode           <- opt[["mode"]]
gene_map_file  <- opt[["gene-map"]]

# Entrez gene ID -> HGNC symbol map (optional; falls back to raw ID when absent)
entrez2symbol <- character(0)
if (!is.null(gene_map_file) && file.exists(gene_map_file)) {
  gm <- read.csv(gene_map_file, colClasses = "character", check.names = FALSE)
  if (all(c("entrez", "symbol") %in% colnames(gm))) {
    gm <- gm[!is.na(gm$symbol) & nchar(gm$symbol) > 0, ]
    entrez2symbol <- setNames(gm$symbol, gm$entrez)
  }
}
to_symbol <- function(ids) {
  ids <- as.character(ids)
  sym <- entrez2symbol[ids]
  ifelse(is.na(sym) | nchar(sym) == 0, ids, sym)
}

# ---------------------------------------------------------------------------
# Publication theme
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
# Load reference UMAP. Prefer the PREDICT-frame reference (UMAP_predict_ref_<mode>.csv,
# written by moma_project.R) so the control cluster + center cross live in the SAME frame
# as the KO dots and WT cross; fall back to the training embedding (UMAP_<mode>_raw_data)
# for legacy results.
# ---------------------------------------------------------------------------
ref_pf_csv <- file.path(output_viz_dir, sprintf("UMAP_predict_ref_%s.csv", mode))
ref_csv    <- file.path(output_viz_dir, sprintf("UMAP_%s_raw_data.csv", mode))
if (file.exists(ref_pf_csv)) {
  ref_csv <- ref_pf_csv
} else if (!file.exists(ref_csv)) {
  stop("Reference UMAP file not found: ", ref_csv)
}

umap_ref <- read.csv(ref_csv, row.names = 1, check.names = FALSE)
ctrl_df  <- umap_ref[umap_ref$condition == "Control", ]
center_x <- median(ctrl_df$UMAP1)
center_y <- median(ctrl_df$UMAP2)

# ---------------------------------------------------------------------------
# Find MOMA result files
# ---------------------------------------------------------------------------
result_files <- list.files(output_viz_dir,
                           pattern   = "^UMAP_Results_.*\\.csv$",
                           full.names = TRUE)
message(sprintf("Found %d MOMA result file(s)", length(result_files)))

if (length(result_files) == 0) {
  message("No MOMA result files found — nothing to plot.")
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Plot each sample
# ---------------------------------------------------------------------------
target_rows <- list()

for (result_file in result_files) {

  bname       <- sub("\\.csv$", "", basename(result_file))
  bname       <- sub("^UMAP_Results_", "", bname)
  sample_name <- sub("^MOMA_target_results_", "", bname)

  tmp_df  <- read.csv(result_file, row.names = 1, check.names = FALSE)
  if (!"recovery_proj_pca" %in% colnames(tmp_df)) {
    stop("recovery_proj_pca column missing in ", basename(result_file),
         " — re-run moma_project.R to regenerate UMAP_Results_*.csv.")
  }
  # Rank targets by recovery toward Control
  # (recovery_proj_pca descending; > 0 = toward Control, larger = stronger).
  top_idx <- order(tmp_df$recovery_proj_pca, decreasing = TRUE)[seq_len(min(10, nrow(tmp_df)))]
  top_df      <- tmp_df[top_idx, ]
  top_df$gene <- to_symbol(rownames(top_df))   # Entrez ID → HGNC symbol

  # Cross marker data frame (Control Center + optional Sample position)
  cross_df <- data.frame(
    x     = center_x,
    y     = center_y,
    group = "Control Center",
    stringsAsFactors = FALSE
  )

  # Sample WT position = the projected WT baseline (wt_x/wt_y from moma_project.R) that
  # the KO dots and the recovery score are anchored to — NOT umap_ref[sample], which is
  # the training-fit position and lands in a different frame (making the cross look
  # inconsistent with the dots). Fall back to umap_ref only for legacy results without
  # wt_x/wt_y.
  sample_label <- sprintf("Sample(%s)", sample_name)
  if (all(c("wt_x", "wt_y") %in% colnames(tmp_df)) && !is.na(tmp_df$wt_x[1])) {
    sx <- tmp_df$wt_x[1]; sy <- tmp_df$wt_y[1]
  } else if (sample_name %in% rownames(umap_ref)) {
    sx <- umap_ref[sample_name, "UMAP1"]; sy <- umap_ref[sample_name, "UMAP2"]
  } else {
    sx <- NA_real_; sy <- NA_real_
  }
  has_sample <- !is.na(sx)
  if (has_sample) {
    cross_df <- rbind(cross_df, data.frame(
      x = sx, y = sy, group = sample_label, stringsAsFactors = FALSE))
  }

  # Colorblind-safe color mapping (Wong 2011)
  color_map <- c(
    "Control"        = "#009E73",   # bluish green — Control reference (Wong)
    "Top-10 targets" = "#E69F00",   # orange — MOMA top targets
    "Control Center" = "#333333"    # dark   — center cross marker
  )
  if (has_sample) {
    color_map[sample_label] <- "#D55E00"   # vermilion — sample WT (before-KO) position
  }

  # Arrow data: sample WT (before-KO) position → KO position for top-10 targets
  arrow_df <- NULL
  if (has_sample) {
    arrow_df <- data.frame(
      x      = sx,
      y      = sy,
      xend   = top_df$X,
      yend   = top_df$Y,
      toward = top_df$recovery_proj_pca > 0
    )
  }

  p <- ggplot() +
    stat_ellipse(data = ctrl_df,
                 aes(x = UMAP1, y = UMAP2, color = "Control", fill = "Control"),
                 type = "norm", level = 0.80, linetype = "dashed", linewidth = 0.3,
                 geom = "polygon", alpha = 0.06, show.legend = FALSE) +
    stat_ellipse(data = ctrl_df,
                 aes(x = UMAP1, y = UMAP2, color = "Control"),
                 type = "norm", level = 0.80, linetype = "dashed", linewidth = 0.3,
                 show.legend = FALSE) +
    geom_point(data = ctrl_df,
               aes(x = UMAP1, y = UMAP2, color = "Control"),
               alpha = 0.7, size = 2.4, shape = 16) +
    # All MOMA gene results (gray, no legend)
    geom_point(data = tmp_df,
               aes(x = X, y = Y),
               color = "gray70", alpha = 0.3, size = 1.4, shape = 16,
               show.legend = FALSE) +
    # Movement arrows: sample → KO position for top-10
    #   vermilion = toward Control (recovery axis), grey = away
    { if (!is.null(arrow_df))
        geom_segment(data = arrow_df[!arrow_df$toward, , drop = FALSE],
                     aes(x = x, y = y, xend = xend, yend = yend),
                     color = "grey70", alpha = 0.40, linewidth = 0.3,
                     arrow = arrow(length = unit(0.12, "cm"), type = "open"),
                     show.legend = FALSE)
      else list() } +
    { if (!is.null(arrow_df))
        geom_segment(data = arrow_df[arrow_df$toward, , drop = FALSE],
                     aes(x = x, y = y, xend = xend, yend = yend),
                     color = "#D55E00", alpha = 0.50, linewidth = 0.3,
                     arrow = arrow(length = unit(0.12, "cm"), type = "open"),
                     show.legend = FALSE)
      else list() } +
    geom_point(data = top_df,
               aes(x = X, y = Y, color = "Top-10 targets"),
               alpha = 0.9, size = 2.4, shape = 16) +
    # Cross markers: Control Center + Sample (shape = 3 is "+")
    geom_point(data = cross_df,
               aes(x = x, y = y, color = group),
               shape = 3, size = 3.5, stroke = 0.9) +
    # Gene labels for top-10 (HGNC symbols, italic; halo + leader lines for legibility)
    geom_text_repel(data = top_df,
                    aes(x = X, y = Y, label = gene),
                    size = 3.0, fontface = "italic", seed = 42,
                    max.overlaps = Inf, force = 8, force_pull = 0.1,
                    box.padding = 0.9, point.padding = 0.5,
                    min.segment.length = 0, segment.size = 0.3,
                    segment.color = "grey55", segment.alpha = 0.8,
                    max.time = 5, max.iter = 200000,
                    bg.color = "white", bg.r = 0.12) +
    scale_color_manual(values = color_map, name = NULL) +
    scale_fill_manual(values = color_map, guide = "none") +
    labs(x = "UMAP 1", y = "UMAP 2",
         title = "MOMA Simulation results") +
    theme_pub() +
    theme(aspect.ratio = 1)   # physically square UMAP panel

  out_png <- file.path(output_viz_dir, sprintf("UMAP_%s.png", bname))
  out_svg <- file.path(output_viz_dir, sprintf("UMAP_%s.svg", bname))
  ggsave(out_png, p, width = 7.2, height = 6.2, dpi = 300, bg = "white")
  ggsave(out_svg, p, width = 7.2, height = 6.2, bg = "white")
  message(sprintf("  Saved: %s", basename(out_png)))

  target_rows[[bname]] <- data.frame(
    Gene               = rownames(top_df),
    Symbol             = top_df$gene,
    Sample             = sample_name,
    rank               = seq_len(nrow(top_df)),
    recovery_proj_pca  = top_df$recovery_proj_pca,
    delta_distance     = top_df$delta_distance,
    distance_to_center = top_df$distance_to_center,
    stringsAsFactors   = FALSE
  )
}

# ---------------------------------------------------------------------------
# Save target_genes.csv
# ---------------------------------------------------------------------------
if (length(target_rows) > 0) {
  target_genes_df <- do.call(rbind, target_rows)
  write.csv(target_genes_df,
            file.path(output_viz_dir, "target_genes.csv"),
            row.names = FALSE)
  message("Saved: target_genes.csv")
}

message("MOMA visualization done.")
