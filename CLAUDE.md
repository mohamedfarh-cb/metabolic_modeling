# tf_metabolism — Project Memory (CLAUDE.md)

Guidance for any AI assistant (Claude Code / Cowork) or human working on this
repository. Read this first. It records what the pipeline does, the local
extensions added on top of the upstream framework, and the non-obvious gotchas
that are easy to get wrong.

---

## 1. What this project is

`tf_metabolism` characterizes **metabolic reprogramming between two conditions**
by integrating transcriptomics with a genome-scale metabolic model (GEM). Given
two transcriptome groups (e.g. control vs disease), it:

1. Finds differentially expressed transcripts (DEGs).
2. Reconstructs condition-specific GEMs from **Recon2M.2** via tINIT.
3. Predicts per-sample metabolic fluxes (Gurobi LP).
4. Computes **differential flux (DE flux)** between conditions.
5. Runs flux-sum and pathway enrichment.
6. Runs **MOMA single-gene knockout** simulation per patient.
7. **(Local extension)** Aggregates MOMA results into a **cohort-level gene
   target ranking**.
8. Produces UMAP/PCA and figure outputs.

Entry point: `tf_metabolism/__main__.py`. Wrapper: `tf_metabolism.py`.

Base command:
```bash
python ./tf_metabolism.py -o <out_dir> -c1 <control.csv> -c2 <test.csv>
```
Inputs `-c1/-c2` are CSVs with Entrez gene IDs as rows, samples as columns.

Environment: Python 3.6, R 4.2, Gurobi 9.1 (license required). See `environment.yml`.

---

## 2. Current biological project context (IBDMDB)

The pipeline is being applied to the IBDMDB (HMP2) cohort, comparing IBD vs
nonIBD with tissue matching:

- **CD-ileum**: 20 nonIBD vs 36 CD (condition2 = CD, n=36).
- **UC-rectum**: 22 nonIBD vs 20 UC (condition2 = UC, n=20).

Inputs are VST-normalized, Entrez-mapped, tissue-matched count matrices
(`*_rectum_vst_entrez.csv`, etc.) from the IBDMDB project, NOT the bundled
example data.

---

## 3. Local extensions added on top of upstream tf_metabolism

These were added in this project and are NOT in the original upstream repo.

### 3a. DE-flux filtering CLI flags (implemented)
The DE-flux thresholds, previously hard-coded, are now CLI flags:
- `--flux_pval`   (default 0.05)
- `--flux_filter` (`raw_p`/`FDR`/`Bonferroni`/`None`, default `raw_p`)
- `--flux_rbc`    (default 0.5) — effect-size gate, |rank-biserial| cutoff

### 3b. Cohort target-ranking module (implemented)
`tf_metabolism/statistical_analysis/cohort_target_ranking.py` aggregates the
per-patient MOMA knockout scores into a cohort-level ranking. Flags:
- `--target_fdr`        (default 0.05)
- `--target_effect_min` (default 0.1)
- `--target_topk`       (default 10)
Outputs `cohort_target_ranking.csv` and `viz/moma_gene_scores_full.csv`.

### 3c. Per-patient score normalization (implemented)
The MOMA recovery score is normalized per patient inside the ranking step
(fixed a bug where raw scores were ~1e9–1e11 and saturated the statistics).

### NOT yet implemented (see FUTURE_REQUESTS.md)
- `--target_method`, `--target_switch_n` flags (size-adaptive method control).
- Option B: effect size on BOTH test branches (Cohen's d for parametric).
- Investigation of why the mixed-effects model never runs.

---

## 4. Key data-flow facts (easy to get wrong)

- **DE flux uses `raw_p` + effect-size gate, NOT FDR.** Reactions enter
  `Differential_fluxes.csv` if raw p < 0.05 AND |rank-biserial| > `--flux_rbc`.
  Nothing in the IBD flux data survives FDR, which is why raw_p + effect size is
  used. **Treat DE flux and everything downstream as exploratory /
  hypothesis-generating, not FDR-significant.**

- **DEG (transcriptome) DE is different from DE flux.** The transcript comparison
  (`Differentially_expressed_transcripts.csv`) uses the function default
  `filter_by='FDR'` and NO effect-size gate. It is unaffected by `--flux_rbc`.

- **`cohort_target_ranking.csv` is built from DE flux, NOT all flux.** The MOMA
  recovery scores are measured in the PCA/UMAP space defined by the
  **differential-flux reactions** (`run_r_umap(..., mode='df')`,
  `export_moma_gene_scores_full(..., mode='df')`). Therefore the ranking is
  **cutoff-dependent**: changing `--flux_rbc` changes the DF reaction set →
  changes the scoring space → changes the ranking. `All_flux_comparison.csv` is a
  separate unfiltered side output and is NOT used by the ranking.

- **MOMA knocks out all model genes**, but each gene's *score* is its recovery
  measured only along the DE-flux reactions. So "candidate genes" = whole model;
  "how good each is" = judged through the DE-flux lens only.

- **The effect-size gate is currently coupled to the test branch.** In
  `statistical_comparison.py`, rank-biserial is computed only on the
  non-parametric (Wilcoxon) path (n < 30). At n ≥ 30 the code switches to the
  t-test and the rbc gate silently disappears. For IBD cohorts (n ≤ 36) this
  mostly stays on the Wilcoxon path, but it is a latent large-n bug — see
  FUTURE_REQUESTS.md "Option B".

---

## 5. Reading `cohort_target_ranking.csv`

One row per gene (GPR-coupled genes share identical scores; `gpr_group` flags
them). A strong target = low `rra_fdr`, high `recurrence_prop`, high
`bootstrap_freq`, `significant = True`.

| Column | Meaning |
|---|---|
| `rra_score` / `rra_fdr` | Robust Rank Aggregation; lower = more consistent target |
| `recurrence_count`/`prop` | How many / fraction of patients rank it top-k |
| `wilcoxon_q` / `rank_biserial` / `median_score` | Effect-consistency test |
| `bootstrap_freq` | Robustness to cohort resampling |
| `perm_pvalue` | Permutation null sensitivity |
| `kendall_W` | Cohort-level ranking agreement (same in every row) |
| `mixedlm_estimate` | Mixed-effects estimate — CURRENTLY ALWAYS NaN (never runs) |
| `significant` | FDR AND effect-size pass |

---

## 6. Resuming a run cheaply (skip the expensive steps)

GEM reconstruction, flux prediction, and MOMA are slow (Gurobi). Each has a
skip-guard. To re-run only DE flux → ranking in a fresh `-o` folder, that folder
must contain:

- `condition1/reconstructed_models/` + `condition1/metabolic_gene_expression.csv`
- `condition2/reconstructed_models/` + `condition2/metabolic_gene_expression.csv`
- `flux1.csv`, `flux2.csv`
- `targeting_results/` (to skip MOMA; cutoff-independent, safe to reuse)

Everything else (`Differential_fluxes.csv`, pathway tables, `viz/`,
`cohort_target_ranking.csv`) is regenerated. Run into a SEPARATE output folder
(e.g. `output_CD_ileum_rbc_0.3`) to avoid overwriting a previous cutoff's results.

**Outputs are gitignored (`output_*/`) and NOT recoverable from git.** Copy aside
before deleting; `rm` in a terminal is permanent (no Trash).

---

## 7. Conventions

- Patient IDs: letter + 4 digits (e.g. C3002, H4004, M2039).
- Output folder naming: `output_<cohort>_<tissue>_<cutoff>` (e.g.
  `output_CD_ileum_rbc_0.3`).
- All randomness seeded (seed 42 in the R/Python dimred and ranking).
- Solver: Gurobi > CPLEX > GLPK.

---

## 8. Pointers

- `README.md` — upstream usage and output descriptions.
- `WORKLOG.md` — session-by-session progress and findings.
- `FUTURE_REQUESTS.md` — known issues and planned improvements.
