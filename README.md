# tf_metabolism

**Characterizing metabolic reprogramming between two conditions**

A framework that characterizes metabolic reprogramming between two conditions by integrating
transcriptomic data with a genome-scale metabolic model (GEM). It reconstructs condition-specific
models, predicts metabolic fluxes, identifies differentially active reactions and the metabolic
pathways that are reprogrammed, and ranks candidate metabolic intervention targets via MOMA
single-gene knockout simulation.

## Overview

The pipeline (`tf_metabolism/__main__.py`) runs end to end:

1. **Differential expression** — two-group comparison of input transcriptomes.
2. **GEM reconstruction** — condition-specific models built from Recon2M.2 via tINIT
   (omics integration in `tf_metabolism/omics_integration/`).
3. **Flux prediction** — context-specific fluxes per sample (`metabolic_simulation/`,
   Gurobi-based LP/MOMA).
4. **Differential flux & flux-sum analysis** — statistically compares fluxes and
   metabolite turnover between conditions (`statistical_analysis/`).
5. **Pathway enrichment** — enrichment analysis over the differentially active metabolic
   reactions to highlight reprogrammed pathways (`statistical_analysis/enrichment.py`).
6. **MOMA targeting simulation** — single gene knockouts on each condition-specific GEM,
   scored by flux perturbation via MOMA (`run_targeting_simulation`, writes `targeting_results/`).
7. **Visualization** — UMAP embeddings and publication figures via the bundled R scripts
   (`flux_dimred.R`, `moma_project.R`, `moma_viz.R`). The enriched-pathway
   bar plot has a selectable backend: R/ggplot2 or matplotlib
   (`visualize_pathway_enrichment(..., use_r=True|False)`, both write `viz/pathway_barplot.png`).
8. **Cross-cohort target ranking** — aggregates the per-patient MOMA rankings into one
   cohort-level ranking of candidate target genes
   (`statistical_analysis/cohort_target_ranking.py`). It is *size-adaptive*, mirroring the
   pipeline's n=30 parametric/non-parametric switch: Robust Rank Aggregation (RRA, pure
   Python) + a within-patient permutation null are primary at all sizes, and a linear
   mixed-effects model is added at n≥30 (graceful RRA-only fallback otherwise). Genes are
   flagged significant only when they pass BOTH an FDR threshold AND a minimum effect size,
   which keeps results meaningful at large n. GPR-collinear genes (identical MOMA scores,
   e.g. PYCR1/PYCR2) are grouped via a `gpr_group` column. Writes
   `cohort_target_ranking.csv` (full per-patient scores in `viz/moma_gene_scores_full.csv`).

## Requirements

- **Python 3.6** and **R 4.2** (the pipeline calls R scripts via `Rscript`)
- **Gurobi 9.1** with a valid license — free academic licenses at
  <https://www.gurobi.com/academia/>
- Developed and tested on Linux (x86-64).

## Installation

### Option A — Conda (recommended, one command)

Installs Python packages, R, R packages, and Gurobi together:

```bash
conda env create -f environment.yml
conda activate gems
```

> Gurobi still requires you to activate a license (`grbgetkey <your-key>`).

### Option B — pip + R (manual)

```bash
# Python
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

# R packages (see the commented list in requirements.txt)
Rscript -e 'install.packages(c("cluster","dplyr","ggplot2","ggrepel","jsonlite",
  "optparse","patchwork","scales","stringr","svglite","umap"))'
```

Make sure `Rscript` is on your `PATH` (or set the `RSCRIPT` environment variable).

## Usage

```bash
python ./tf_metabolism.py -o <output_dir> -c1 <condition1.csv> -c2 <condition2.csv>
```

Example (the configuration in `run.sh`, healthy control vs. Crohn's disease):

```bash
python ./tf_metabolism.py -o ./output_CD \
    -c1 ./input_data/HC_entrez_tpm.csv \
    -c2 ./input_data/CD_entrez_tpm.csv
```

### Arguments

| Flag | Argument | Description |
|------|----------|-------------|
| `-o`  / `--output_dir`        | path | Output directory (created if missing) |
| `-c1` / `--input_omics_file1` | CSV  | Transcriptome for condition 1 (genes × samples) |
| `-c2` / `--input_omics_file2` | CSV  | Transcriptome for condition 2 (genes × samples) |
| `--flux_pval`   | float | DE-flux p-value cutoff (default `0.05`) |
| `--flux_filter` | `raw_p` \| `FDR` \| `Bonferroni` \| `None` | DE-flux p-value filter mode (default `raw_p`) |
| `--flux_rbc`    | float | DE-flux effect-size cutoff \|rank-biserial\| (default `0.5`) |
| `--target_fdr`        | float | Cross-cohort target-ranking BH-FDR cutoff (default `0.05`) |
| `--target_effect_min` | float | Min effect size (recurrence proportion) for a significant target (default `0.1`) |
| `--target_topk`       | int   | Top-k cutoff for recurrence / bootstrap selection (default `10`) |
| `--target_method`     | `auto` \| `rra` \| `mixedlm` | Target-ranking inference: `auto` size-adaptive, `rra` RRA only, `mixedlm` force mixed model (default `auto`) |
| `--target_switch_n`   | int   | Patient count at which `auto` adds the mixed-effects model (default `40`) |
| `--target_nperm`      | int   | Permutation-null resamples for target ranking (default `1000`) |
| `--target_nboot`      | int   | Bootstrap resamples for target ranking (default `1000`) |
| `--ko_targets`        | GENE [GENE ...] | Render a cohort-level knockout figure per gene (`viz/ko_cohort_<SYMBOL>.png` + `.svg`): every condition2 patient's wild-type → knockout shift against the nonIBD control cloud. Accepts HGNC symbols or model gene IDs, space- or comma-separated. Omit to skip. Requires a completed MOMA projection step; renders from existing outputs without re-running MOMA. Display only — MOMA targets are hypothesis-only (see FUTURE_REQUESTS [P0]) (default: none) |

The `--flux_*` flags tune the differential-flux filtering in step 4
(`Differential_fluxes.csv`). The defaults reproduce the original behavior, so
omitting them leaves output unchanged. For example, to apply FDR correction and a
stricter effect-size threshold:

```bash
python ./tf_metabolism.py -o ./output_CD \
    -c1 ./input_data/HC_entrez_tpm.csv \
    -c2 ./input_data/CD_entrez_tpm.csv \
    --flux_pval 0.01 --flux_filter FDR --flux_rbc 0.7
```

The `--target_*` flags control the cross-cohort target ranking in step 8
(`cohort_target_ranking.csv`). RRA + permutation are the always-on backbone in
every mode; `--target_method` only governs the mixed-effects model — `auto`
(default) adds it once the cohort reaches `--target_switch_n` patients (and it
converges), `rra` never fits it, and `mixedlm` forces it. The defaults reproduce
the original behavior. For example, to force RRA-only ranking with heavier
resampling:

```bash
python ./tf_metabolism.py -o ./output_CD \
    -c1 ./input_data/HC_entrez_tpm.csv \
    -c2 ./input_data/CD_entrez_tpm.csv \
    --target_method rra --target_nperm 5000 --target_nboot 5000
```

### Input format

- `-c1` / `-c2`: CSV with gene identifiers (Entrez) as the row index and samples as columns.
  See `input_data/HC_entrez_tpm.csv`, `input_data/CD_entrez_tpm.csv`,
  or the smaller `input_data/sample1.csv` / `sample2.csv` examples.

### Outputs

Written under `<output_dir>/`:

| File / directory | Contents |
|------------------|----------|
| `condition1/`, `condition2/`   | Reconstructed condition-specific GEMs and split omics |
| `flux1.csv`, `flux2.csv`       | Predicted flux profiles per condition |
| `Differentially_expressed_transcripts.csv` | DEGs between the two conditions |
| `Differential_fluxes.csv`, `All_flux_comparison.csv` | Differential metabolic fluxes |
| `flux_sum_results.csv`         | Flux-sum (metabolite turnover) comparison |
| `gene_subsystem_flux.csv`      | Gene-to-subsystem flux mapping |
| `Up_regulated_pathways.csv`, `Down_regulated_pathways.csv` | Enriched metabolic pathways |
| `targeting_results/`           | MOMA single-gene knockout flux perturbations |
| `cohort_target_ranking.csv`    | Cohort-level ranking of candidate target genes (size-adaptive stats) |
| `viz/`                         | UMAP embeddings, `pathway_barplot.png`, `target_oncoprint.png` (+ .svg), `moma_gene_scores_full.csv`, figures |

## Repository layout

```
tf_metabolism/
├── __main__.py                # pipeline entry point
├── metabolic_model/           # GEM editing, GPR manipulation
├── metabolic_simulation/      # flux prediction, MOMA, gap filling, tasks
├── omics_integration/         # tINIT / INIT model reconstruction
├── statistical_analysis/      # differential comparison, enrichment, cohort target ranking
├── utils.py                   # R-script orchestration, UMAP helpers
└── r_based_viz/               # R visualization scripts (ggplot2 / umap)
    ├── flux_dimred.R           #   PCA + UMAP embedding
    ├── moma_project.R          #   project MOMA results into UMAP/PCA space
    ├── moma_viz.R              #   MOMA UMAP figures + target ranking
    └── pathway_enrichment_plot.R  # enriched-pathway bar plot
input_data/                    # example transcriptomes, model, media, tasks
environment.yml                # conda environment (Python + R + Gurobi)
requirements.txt               # Python deps (+ R deps as comments)
```
