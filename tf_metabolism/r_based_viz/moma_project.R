#!/usr/bin/env Rscript
# =============================================================================
# moma_project.R — Project MOMA simulation results into R UMAP space
#
# Uses PCA + UMAP models saved by flux_dimred.R (RDS files).
# Delta-based correction:
#   projected = merged_flux[patient] + (moma_flux - raw_flux2[patient])
#
# Output columns (per gene KO/KD):
#   X, Y                  — projected UMAP coords
#   distance_to_center    — Euclidean dist to Control center (UMAP, legacy)
#   wt_distance_to_center — sample's own WT dist to Control center (reference)
#   delta_distance        — d(KO) - d(WT);   < 0 = moved TOWARD Control
#   recovery_proj_umap    — displacement·(Test->Control axis) in UMAP; > 0 = TOWARD
#   recovery_proj_pca     — same in PCA space (linear, preferred for scoring)
#
# Usage:
#   Rscript moma_project.R \
#       --output-dir     <path>    # contains flux2.csv, targeting_results/
#       --output-viz-dir <path>    # contains *_model.rds, merged_flux.csv
#       [--mode          df]       # default: df
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(umap)
  library(jsonlite)
})

option_list <- list(
  make_option("--output-dir",     type = "character",
              help = "dir containing flux2.csv and targeting_results/"),
  make_option("--output-viz-dir", type = "character",
              help = "dir containing R model RDS files and merged_flux.csv"),
  make_option("--mode",           type = "character", default = "df",
              help = "UMAP mode (default: df)"),
  make_option("--moma-subdir",    type = "character", default = "",
              help = "subdir under targeting_results/ to read MOMA files from, and under output-viz-dir/ to write results to (default: empty = root)"),
  make_option("--file-prefix",    type = "character", default = "MOMA_target_results_",
              help = "filename prefix for MOMA result CSVs (default: MOMA_target_results_)")
)

opt            <- parse_args(OptionParser(option_list = option_list))
output_dir     <- opt[["output-dir"]]
output_viz_dir <- opt[["output-viz-dir"]]
mode           <- opt[["mode"]]
moma_subdir    <- opt[["moma-subdir"]]
file_prefix    <- opt[["file-prefix"]]

# ---------------------------------------------------------------------------
# Load R models
# ---------------------------------------------------------------------------
pca_rds  <- file.path(output_viz_dir, sprintf("pca_%s_model.rds",  mode))
umap_rds <- file.path(output_viz_dir, sprintf("umap_%s_model.rds", mode))
if (!file.exists(pca_rds))  stop("PCA model not found: ",  pca_rds)
if (!file.exists(umap_rds)) stop("UMAP model not found: ", umap_rds)

pca_model  <- readRDS(pca_rds)
umap_model <- readRDS(umap_rds)

params_json <- file.path(output_viz_dir, sprintf("umap_%s_params.json", mode))
params      <- fromJSON(params_json)
n_pcs       <- as.integer(params$n_pcs)
features    <- params$features

# ---------------------------------------------------------------------------
# Reference: Control center in R UMAP space
# ---------------------------------------------------------------------------
umap_ref <- read.csv(file.path(output_viz_dir, sprintf("UMAP_%s_raw_data.csv", mode)),
                     row.names = 1, check.names = FALSE)
ctrl_df  <- umap_ref[umap_ref$condition == "Control", ]
center_x <- median(ctrl_df$UMAP1)
center_y <- median(ctrl_df$UMAP2)

# Test center + recovery axis in UMAP space (Test -> Control), unit vector
test_df_ref  <- umap_ref[umap_ref$condition == "Test", ]
test_x       <- median(test_df_ref$UMAP1)
test_y       <- median(test_df_ref$UMAP2)
rec_umap_vec <- c(center_x - test_x, center_y - test_y)
rec_umap_vec <- rec_umap_vec / sqrt(sum(rec_umap_vec^2))

# ---------------------------------------------------------------------------
# Load merged flux baseline and raw wild-type flux
# ---------------------------------------------------------------------------
merged_df <- read.csv(file.path(output_viz_dir, "merged_flux.csv"),
                      row.names = 1, check.names = FALSE)

flux2_df  <- read.csv(file.path(output_dir, "flux2.csv"),
                      row.names = 1, check.names = FALSE)
flux2_df[is.na(flux2_df)] <- 0

# ---------------------------------------------------------------------------
# PCA transform helper (applies training center/scale then rotates)
# ---------------------------------------------------------------------------
project_pca <- function(X_mat, pca_model, n_pcs) {
  # X_mat: n_genes × n_features  (rows = new observations, cols = features)
  X_c <- sweep(X_mat, 2, pca_model$center, "-")
  if (!is.null(pca_model$scale) && !isFALSE(pca_model$scale)) {
    X_c <- sweep(X_c, 2, pca_model$scale, "/")
  }
  X_c %*% pca_model$rotation[, seq_len(n_pcs), drop = FALSE]
}

# ---------------------------------------------------------------------------
# Recovery axis in PCA space (Test -> Control centroid), unit vector
#   Linear & faithful to flux distances — preferred over 2D UMAP for scoring.
#   Centroids = median of each group's ComBat baseline projected into PC space.
# ---------------------------------------------------------------------------
ctrl_ids <- intersect(rownames(umap_ref)[umap_ref$condition == "Control"],
                      colnames(merged_df))
test_ids <- intersect(rownames(umap_ref)[umap_ref$condition == "Test"],
                      colnames(merged_df))
ctrl_base <- t(merged_df[features, ctrl_ids, drop = FALSE])   # samples × features
test_base <- t(merged_df[features, test_ids, drop = FALSE])
ctrl_pca  <- project_pca(ctrl_base, pca_model, n_pcs)          # samples × n_pcs
test_pca  <- project_pca(test_base, pca_model, n_pcs)
ctrl_centroid_pca <- apply(ctrl_pca, 2, median)
test_centroid_pca <- apply(test_pca, 2, median)
rec_pca_vec <- ctrl_centroid_pca - test_centroid_pca
rec_pca_vec <- rec_pca_vec / sqrt(sum(rec_pca_vec^2))

# ---------------------------------------------------------------------------
# Project each MOMA result file
# ---------------------------------------------------------------------------
moma_read_dir <- if (nchar(moma_subdir) > 0)
  file.path(output_dir, "targeting_results", moma_subdir) else
  file.path(output_dir, "targeting_results")

umap_out_dir <- if (nchar(moma_subdir) > 0)
  file.path(output_viz_dir, moma_subdir) else output_viz_dir
dir.create(umap_out_dir, showWarnings = FALSE, recursive = TRUE)

file_pattern <- sprintf("^%s.*\\.csv$", file_prefix)
moma_files <- list.files(moma_read_dir,
                         pattern   = file_pattern,
                         full.names = TRUE)
message(sprintf("[MOMA-Project-%s] Projecting %d files from %s ...",
                mode, length(moma_files), moma_read_dir))

for (moma_file in moma_files) {
  base_noext  <- sub("\\.csv$", "", basename(moma_file))
  sample_name <- sub(sprintf("^%s", file_prefix), "", base_noext)

  df_moma <- read.csv(moma_file, row.names = 1, check.names = FALSE)  # reactions × genes

  # Reindex rows to features (fill missing with 0), transpose to genes × features
  feat_order <- features
  df_feat    <- df_moma[match(feat_order, rownames(df_moma)), , drop = FALSE]
  rownames(df_feat) <- feat_order
  df_feat[is.na(df_feat)] <- 0
  moma_mat <- t(df_feat)   # genes × features

  # Delta correction
  has_wt <- sample_name %in% colnames(merged_df) &&
            sample_name %in% colnames(flux2_df)
  if (has_wt) {
    wt_combat <- as.numeric(merged_df[feat_order, sample_name])
    wt_raw    <- as.numeric(flux2_df[feat_order, sample_name])
    wt_raw[is.na(wt_raw)] <- 0
    delta   <- sweep(moma_mat, 2, wt_raw,    "-")
    X_input <- sweep(delta,    2, wt_combat, "+")
  } else {
    message(sprintf("  WARNING: %s not in combat/flux2 — raw projection", sample_name))
    X_input <- moma_mat
  }

  # PCA → UMAP (perturbed KO/KD points)
  pca_scores  <- project_pca(X_input, pca_model, n_pcs)   # genes × n_pcs
  umap_coords <- predict(umap_model, pca_scores)            # genes × 2

  d_ko        <- sqrt((umap_coords[, 1] - center_x)^2 +
                      (umap_coords[, 2] - center_y)^2)
  gene_names  <- colnames(df_moma)

  if (has_wt) {
    # WT reference projected through the SAME pipeline (delta = 0) for a
    # self-consistent displacement vector (KO and WT both via predict()).
    wt_mat  <- matrix(wt_combat, nrow = 1, dimnames = list("WT", feat_order))
    wt_pca  <- project_pca(wt_mat, pca_model, n_pcs)        # 1 × n_pcs
    wt_umap <- as.numeric(predict(umap_model, wt_pca))      # length 2
    d_wt    <- sqrt((wt_umap[1] - center_x)^2 + (wt_umap[2] - center_y)^2)

    # Recovery-axis projections (> 0 = displaced TOWARD Control)
    disp_umap     <- sweep(umap_coords, 2, wt_umap,            "-")  # genes × 2
    rec_proj_umap <- as.numeric(disp_umap %*% rec_umap_vec)
    disp_pca      <- sweep(pca_scores,  2, as.numeric(wt_pca), "-")  # genes × n_pcs
    rec_proj_pca  <- as.numeric(disp_pca %*% rec_pca_vec)

    result_df <- data.frame(
      X = umap_coords[, 1],
      Y = umap_coords[, 2],
      distance_to_center    = d_ko,
      wt_distance_to_center = d_wt,
      wt_x = wt_umap[1],                      # projected WT baseline coords — the SAME
      wt_y = wt_umap[2],                      # frame as the KO dots; the sample cross must
                                              # be drawn here (NOT umap_ref[sample], which is
                                              # the training-fit position and disagrees).
      delta_distance        = d_ko - d_wt,    # < 0 = toward Control (per-sample)
      recovery_proj_umap    = rec_proj_umap,  # > 0 = toward Control (UMAP 2D)
      recovery_proj_pca     = rec_proj_pca,   # > 0 = toward Control (PCA, preferred)
      row.names = gene_names
    )
  } else {
    result_df <- data.frame(
      X = umap_coords[, 1],
      Y = umap_coords[, 2],
      distance_to_center    = d_ko,
      wt_distance_to_center = NA_real_,
      wt_x = NA_real_,
      wt_y = NA_real_,
      delta_distance        = NA_real_,
      recovery_proj_umap    = NA_real_,
      recovery_proj_pca     = NA_real_,
      row.names = gene_names
    )
  }

  out_csv <- file.path(umap_out_dir, sprintf("UMAP_Results_%s.csv", base_noext))
  write.csv(result_df, out_csv, row.names = TRUE)
}

message(sprintf("[MOMA-Project-%s] Done.", mode))
