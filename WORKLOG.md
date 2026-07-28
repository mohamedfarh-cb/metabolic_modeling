# WORKLOG — tf_metabolism

Session-by-session log of changes, runs, and findings. Newest at top.

---

## Session: 2026-07-29b (per-target MOMA cohort view — new core pipeline step)

Added `--ko_targets`: one figure per knocked-out gene showing every condition2
patient's WT→KO shift against the control cloud, inverting the per-patient MOMA view.
Design spec: `docs/superpowers/specs/2026-07-29-per-target-moma-cohort-view-design.md`;
plan: `docs/superpowers/plans/2026-07-29-per-target-moma-cohort-view.md`.

- **No re-simulation.** `moma_project.R` already projects EVERY gene for EVERY patient;
  the new step re-slices `viz/UMAP_Results_*.csv` by gene instead of by patient. So it
  renders with no Gurobi re-run for any completed run whose MOMA projection came from the
  current **predict-frame** `moma_project.R` — in practice the **two Human2 runs**
  (`output_{CD_ileum,UC_rectum}_Human2_tINIT_Human2_tINIT`). `moma_project.R` (and the
  [P0b] fix in it) is untouched.
- **NOT the two Recon2M.2 runs.** `output_CD_ileum_rbc_0.3` / `output_UC_rectum_rbc_0.3`
  predate commit `e388856` and have neither `viz/UMAP_predict_ref_df.csv` nor `wt_x`/`wt_y`
  columns, so the step **fails loud** on them (verified: `Predict-frame reference not
  found`, exit 1) — correct behaviour, not a bug. `moma_viz.R` carries a legacy-frame
  fallback for those older results; `moma_cohort_viz.R` deliberately does not, since the
  entire figure is a WT→KO displacement and there is no correct legacy frame for it.
  Rendering there requires re-running the MOMA projection step.
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

---

## Session: 2026-06-24e (display fix — KO plots must fit_transform, not transform)

### Root cause of the KO-plot collapse (it is NOT UMAP vs PCA)
Confirmed in code (tf_metabolism/utils.py):
- `run_r_umap` (L328) builds the baseline UMAP with **fit_transform** (baselines only) — fine.
- `run_moma_results_umap` (L536, per-patient MOMA plots) places KO points with
  **reducer.transform()** — OOD collapse.
- `plot_gene_ko_shift` (L1254, ko_shift plot) also places KO points with
  **reducer.transform()** — OOD collapse.
So BOTH KO-plotting functions project knockout fluxes into a baseline-ONLY embedding
via transform(); knockout states are out-of-distribution → they collapse to a corner
(active genes) or pile at the WT point (silent genes, ~98%). This is gene-independent
(FASN==DBT) → a visualization artifact, not biology.

The reference IBD project's FASN UMAP scatters because it **fit the UMAP on the KO
cloud** (fit_transform, KOs in-distribution). Same algorithm, correct construction.
=> Fix BOTH functions to fit_transform on baselines + KO points together. UMAP is
fine; the bug was transform-on-OOD.

### Important boundary (display vs significance)
- A correctly-built UMAP (fit including KOs) is a fine DISPLAY of the per-patient /
  cohort KO shift.
- It is NOT the significance test. Target calls come from the high-D **competitive
  null** (2026-06-24d), never from 2D proximity. Keep the two separate.

### Action
- Claude Code prompt (2026-06-24e) fixes both KO plots to fit_transform, removes the
  superseded absolute-recovery gating, keeps competitive null as the significance test.

---

## Session: 2026-06-24d (FINAL — competitive null is the reliable test; CD has no target)

### The decisive problem: global drift confounds every absolute test
Almost any knockout nudges most patients slightly toward control (cutting flux
pulls the hyperactive disease state toward the lower-flux control region). So every
ABSOLUTE recovery test is confounded:
- "recovery > 0" → flags ~everything;
- `frac_closer>0.5` / `proj_vs_distance_agree` majority gate → knife-edge at n=36
  (DBT closer in 17/36=0.472 → False; SDSL 19/36=0.528 → True; difference = 2
  patients; binomial vs coin-flip P≈0.43–0.69 for BOTH → neither beats chance);
- mean-recovery ranking → single-patient-outlier driven (top-mean genes help 1/36).

### The reliable test: competitive null (MTA score + Reporter-Metabolite background)
1. Score each KO for movement toward the control centroid (MTA/rMTA-style,
   directional, 30-PC PCA space).
2. Assess significance vs a RANDOM/sampled-knockout background (Reporter-Metabolite
   restandardization: subtract mean / divide SD of random KOs of matched magnitude),
   FDR-corrected. Target = recovers MORE than a random KO → removes the drift.
This is GEM-standard, not ad-hoc: Patil & Nielsen 2005 PNAS (reporter metabolites,
random-set background); Yizhak MTA 2013 / Valcárcel rMTA 2019 (disease→healthy KO
ranking, relative top-10-20%); competitive-vs-self-contained null (gene-set theory).

### Prototype result (drift-controlled permutation null, existing CD data, n=36)
- **0 / 1362 genes pass p<0.05; 0 pass FDR<0.05.**
- DBT, FASN, SDSL, BCKDHA, HADHA, DHFR, FPGS, OLAH — all **p ≈ 0.77** (= random).
- **FINAL CD CONCLUSION: no single-gene knockout robustly reverts the CD metabolic
  state beyond a non-specific perturbation effect.** Supersedes ALL prior shortlists
  (BCAA/FAO from the magnitude score; DHFR/FASN/folate from the raw-distance proto).

### Why DBT etc. looked good on the old metrics (reconciliation)
- Per-patient "closer": real but it is the global drift — DBT is in the majority,
  not above it.
- RRA/recurrence "best": reward consistency of being a LARGE, top-ranked
  perturbation; DBT is near-essential → big consistent effect → high RRA. They
  measure magnitude+consistency, NOT specificity vs background.
- Only the competitive null asks "more than a random KO" → DBT is average (p=0.77).

### Visualization note (figures are NOT evidence of recovery)
- ko_shift / per-patient MOMA UMAPs collapse all "active" KOs to one corner because
  they `umap.transform()` out-of-distribution KO fluxes through a baseline-only
  embedding; ~98% of KOs are silent (pile at the WT point). The collapse is
  gene-independent (FASN == DBT). The reference IBD project's nicer FASN scatter
  used `fit_transform` (KOs in the embedding), not transform — a plotting choice,
  NOT a biology or reliability difference. Do not judge recovery from these 2D plots.

### Action
- Implement competitive null as the official `significant` rule (Claude Code prompt
  written 2026-06-24d). Clean up the superseded score/gate code and stale figures.
- [P0] in FUTURE_REQUESTS.md updated; intermediate attempts marked superseded.

---

## Session: 2026-06-24c (literature check — which recovery metric is conventional)

Question: for MOMA KO simulations, is `delta_distance` or `recovery_proj_pca` the
conventional recovery readout? Findings (web search):
- **MOMA (Segrè, Vitkup & Church 2002, PNAS)** defines only the WT↔KO flux distance
  (min Σ(v_KO−v_WT)²). It has NO "toward control" concept — recovery is an add-on.
- **Signature-reversal field prefers a normalized DIRECTIONAL projection** onto the
  disease→control axis (Signature Projection Score; TReD), explicitly **preferred
  over raw Euclidean distance** because it captures direction AND magnitude relative
  to control.
- **Cancer GEM KO papers** (PNAS 2025 drug-resistant breast cancer: KO + UMAP +
  k-means; redox radiation-resistant tumors) judge recovery by **cluster/region
  membership**, not a raw scalar distance.
Conclusion: `recovery_proj_pca`'s CONCEPT (directional projection) is the standard;
the bug is that this pipeline's version is UNNORMALIZED (~1e11) → magnitude artifact.
Fix = normalized directional projection (primary) + `delta_distance` cross-check.
[P0] in FUTURE_REQUESTS.md updated accordingly.

Refs:
- Segrè et al. 2002, PNAS — https://www.pnas.org/doi/10.1073/pnas.232349399
- Drug-resistant breast cancer KO+UMAP, PNAS 2025 — https://www.pnas.org/doi/10.1073/pnas.2425384122
- Redox radiation-resistant tumor GEM targets — https://www.sciencedirect.com/science/article/pii/S2405471220305020
- Transcriptome signature reversion validation, Brief Bioinform 2023 — https://academic.oup.com/bib/article/24/1/bbac490/6850563
- TReD, Bioinformatics 2025 — https://academic.oup.com/bioinformatics/article/41/10/btaf498/8250709
- Signature Projection Score (direction+magnitude vs Euclidean), USPTO 7,519,519

---

## Session: 2026-06-24b (CRITICAL — MOMA target score measures disruption, not recovery)

### Finding (do NOT carry the BCAA/FAO "targets" forward as validated)
The CD "top targets" (DBT, BCKDHA, BCKDHB, HADHA/B, ECHS1, ACADM) are a SCORING
ARTIFACT, not validated therapeutic targets. The MOMA target score in use ranks
genes by the MAGNITUDE of flux disruption, not by whether the knockout actually
moves CD patients toward the nonIBD (control) state.

Evidence (CD-ileum rbc=0.3, n=36):
- The per-patient score fed to RRA is `recovery_proj_pca` on a raw ~1e11 scale
  (axis projection dominated by displacement magnitude). DBT scores #1 with
  recovery_proj_pca > 0 in 36/36 patients.
- But by the HONEST metric — actual reduction in Euclidean distance to the control
  centroid (`delta_distance` = wt_distance − ko_distance) — DBT barely helps:
  distance_to_center 2.85 → 2.70, and only ~19/36 (53%) of patients move closer.
- The DBT KO-shift UMAP confirms it: all CD points collapse to a tight degenerate
  cluster AWAY from the control cloud (a big, consistent shift, but not toward
  health). DBT is a near-essential central BCAA enzyme → its KO causes a large
  disruption, which the magnitude score rewards.

Re-ranking by TRUE recovery (delta_distance, prototype on existing UMAP_Results):
- DBT/BCKDH drop to rank ~46–50/1362.
- New consistent-recovery top genes: **DHFR, FASN, OLAH, FPGS** (folate/one-carbon
  + de novo lipogenesis), each moving 72% of patients closer; then SDS/SDSL (67%).
  Note FASN matches the user's original example gene.
- Sobering bottom line: NO gene robustly recovers the cohort — max frac_closer =
  0.72; only 4 genes > 70%, 40 genes > 50%. CD's targets are modest/partial, not
  the slam-dunk the magnitude ranking implied.
- Ranking by MEAN recovery is ALSO broken (outlier-driven): the highest-mean genes
  (DCT, DAO, TYRP1) help only 1/36 patients. Use MEDIAN + frac_closer, not mean.

### Fix scope (no expensive rerun needed)
- **Required:** `export_moma_gene_scores_full` (utils.py, the RECENTLY-added
  cohort-ranking Deliverable 1) — change `score` from `recovery_proj_pca` to the
  signed proximity gain `wt_distance_to_center − distance_to_center`. Both columns
  already exist in `viz/UMAP_Results_*.csv` (computed by the ORIGINAL, untouched
  `run_moma_results_umap` / `moma_project.R`). Everything downstream reads
  `moma_gene_scores_full.csv`, so this one change propagates.
- **Optional:** `cohort_target_ranking.py` — add a directional-consistency gate to
  `significant` (require frac_closer > 0.5 and use median signed recovery, not
  mean) to close the single-patient-outlier loophole.
- Not touched: MOMA simulation and the projection layer (already correct). Can
  regenerate `moma_gene_scores_full.csv` + `cohort_target_ranking.csv` from existing
  `UMAP_Results_*.csv` with no Gurobi rerun.

### Action items
- Implement the score fix (Claude Code prompt written 2026-06-24).
- Re-run CD + UC rankings; re-interpret. Expect the BCAA/FAO story to weaken.
- See FUTURE_REQUESTS.md [P0].

---

## Session: 2026-06-24 (gene-ranking statistics — rationale + method wiring)

### What was done
- Wired the remaining cohort target-ranking CLI flags:
  `--target_method {auto,rra,mixedlm}` (default `auto`), `--target_switch_n`
  (default 40), `--target_nperm` (1000), `--target_nboot` (1000). The earlier
  "unrecognized arguments" issue is resolved — all flags now parse and thread
  through to `rank_targets_across_cohort`.
- Decoupled the effect-size gate from the test choice in
  `statistical_comparison.two_grouped_data_comparison` (the "Option B" item):
  the parametric (t-test) branch now reports **Cohen's d** with its own
  `cohens_d_cutoff`, alongside the non-parametric **rank-biserial** /
  `rbc_cutoff`. Generic `Effect size` + `Effect size type` columns record which
  applies. DE-flux call now passes `cohens_d_cutoff=0.8` (the parametric analogue
  of rbc=0.5) so the effect-size filter never silently disappears at n≥30.
- Refined the per-patient score normalization to **median-center + SD-scale**
  (see rationale below); `median_score` is now modest (CD abs-max ≈ 7.5, was
  ~2.3e9) and `rank_biserial` is no longer saturated at ±1 (96% interior).

### Why these statistical methods — gene ranking from MOMA results
The goal is to turn **per-patient** MOMA single-gene-knockout simulations into a
**cohort-level** ranking of candidate targets that is *population/cohort-specific*
and stays valid from small (n<30) to large (n≫100) cohorts. Each method below
answers one failure mode of "just average the per-patient scores."

1. **Per-patient robust scaling (median-center + SD-scale), applied first.**
   Raw MOMA recovery scores are unbounded fluxes (~1e9–1e11) with a large shared
   baseline; a few large-flux reactions dominate, and testing "score vs 0" is
   meaningless because 0 is not the baseline. Centering on the patient **median**
   (the control-centroid baseline) re-expresses each score as *signed recovery
   relative to that baseline* (so the sign/rank tests stop saturating), and the
   **SD** divisor — inflated by the same outliers — keeps even the strongest
   targets on a modest scale (IQR scaling re-inflated them because the bulk genes
   are near-tied, IQR≈44 vs SD≈1.9e10). The transform is **monotone within a
   patient**, so it deliberately leaves all rank-based statistics (RRA,
   recurrence, Kendall's W) unchanged and only fixes the magnitude-based ones
   (`median_score`, mixed-model estimates) and the centering of the signed test.

2. **Robust Rank Aggregation (RRA, Kolde et al. 2012) — primary ranking.**
   Cross-patient flux magnitudes are not comparable even after scaling, but
   *within-patient ranks* are. RRA aggregates each patient's gene ranking via a
   beta order-statistic null, rewarding genes ranked consistently near the top
   across patients and being robust to a minority of patients where a gene looks
   strong by chance. Pure-Python (no new R dependency).

3. **Kendall's W agreement gate.** Before trusting any aggregate, quantify how
   concordant the per-patient rankings actually are. We report the **magnitude**
   of W (not its p-value, which is meaningless at large n) and flag low agreement
   — a cohort with W≈0 has no shared target signal to rank, so a "significant"
   gene there would be an artifact.

4. **Recurrence (top-k count vs Binomial(n, k/G), one-sided).** RRA can rank a
   gene highly on consistently-moderate ranks that are *never actually top-k in
   any patient*. Recurrence is the interpretable, clinically-meaningful effect
   size — "in what fraction of patients is this gene a top target" — and guards
   against that RRA-vs-recurrence disconnect.

5. **Wilcoxon signed-rank + rank-biserial, on the centered signed recovery.**
   Non-parametric (matches the repo's n=30 convention) test of whether a gene's
   recovery is consistently > baseline across patients, with a bounded effect
   size (rank-biserial ∈ [−1,1]). Operates on the median-centered values from (1)
   so it is no longer dominated by the raw-magnitude offset.

6. **Bootstrap over patients (selection frequency).** Reports how often a gene
   stays in the top-k when the patient set is resampled — i.e., robustness to
   *which* patients happened to be sampled, the key concern for a cohort claim.

7. **Within-patient label-permutation null.** Empirical p-values for the
   aggregate statistic under the null "gene labels are exchangeable within a
   patient," controlling for the ranking machinery itself rather than assuming a
   parametric null.

8. **Size-adaptive switch (`method`/`switch_n`).** RRA + permutation are the
   always-on, distribution-free backbone (robust when underpowered). Only for
   large enough cohorts (n ≥ `switch_n`, default 40) does `auto` *add* a
   linear mixed-effects model (per-gene fixed effect, random patient intercept)
   for model-based estimates — with graceful fallback to RRA on non-convergence
   or when the gene count exceeds `mixedlm_max_genes` (300). This mirrors the
   pipeline's existing parametric/non-parametric n=30 philosophy.

9. **Effect-size gate on significance, at every n.** A gene is significant only
   if `rra_fdr < fdr_cutoff` **AND** an effect size (recurrence proportion ≥
   `effect_size_min`) clears a floor. This is the single rule that keeps results
   meaningful at large n, where p-values collapse toward 0 for trivially small
   effects. The decoupled Cohen's-d/rank-biserial work in
   `two_grouped_data_comparison` enforces the same principle upstream in DE-flux.

10. **GPR non-independence grouping.** Genes sharing a reaction (GPR coupling,
    e.g. PYCR1/PYCR2, or the BCKDH subunits DBT/BCKDHA/BCKDHB) get identical MOMA
    scores and would be double-counted as independent evidence. Identical-score
    genes are tagged with a shared `gpr_group` so they are read together.

**Net design intent:** no single statistic is trusted alone. A *cohort-specific*
target must survive concordance (W), consistent ranking (RRA), recurrence,
signed-effect consistency (Wilcoxon), patient-resampling (bootstrap), and a
permutation null — and then clear an explicit effect-size floor. This is why CD
yields a coherent BCAA-catabolism / β-oxidation signal while UC (diffuse,
patient-specific DE flux, n=20) yields none: the methods are built to *reject*
non-reproducible signals, not just to sort them.

### Status of previously-open items (now resolved)
- `--target_method` / `--target_switch_n` (+ `--target_nperm`, `--target_nboot`)
  **implemented** and parsing.
- Effect size on **both** test branches (Cohen's d + rank-biserial) **implemented**.
- Per-patient normalization **refined** (median-center + SD-scale) and verified.

### Still open / by-design notes
- **Mixed-effects model does not fire on the real CD/UC runs** — now for two
  explicit reasons: (a) `switch_n`=40 > CD n=36, and (b) G≈1362 genes >
  `mixedlm_max_genes`=300, so it is skipped with a logged warning and falls back
  to RRA. Both cohorts are therefore scored by the same RRA + permutation
  backbone (as intended), so the CD/UC difference is a real signal difference,
  not a method artifact.
- **RRA-vs-recurrence disconnect** is mitigated by reporting recurrence as the
  effect size and gating on it, but RRA can still surface low-recurrence genes in
  its raw ordering — read `recurrence_prop` alongside `rra_score`.

### Target ranking is cutoff-dependent and NON-monotone (important caveat)
Observed: at rbc=0.5 SLC7A1 was significant=True and DBT was not; at rbc=0.3 the
reverse — DBT (with BCKDHA/BCKDHB) is significant=True (rra_fdr≈3e-54,
recurrence=0.61) and SLC7A1 drops to False. This is expected behaviour, not a bug.

Key distinction (why relaxing the cutoff does NOT simply "add" targets):

- **DE-flux reactions ARE monotone/nested.** A reaction passes if
  `|rank-biserial| > rbc` AND `p<0.05`; since `{|rbc|>0.5} ⊂ {|rbc|>0.3}`, the
  rbc=0.5 reaction set (26 for CD) is a strict subset of the rbc=0.3 set (116).
  Relaxing the cutoff can only add reactions. The user's intuition holds here.

- **Gene targets are NOT monotone.** Targets are not selected by thresholding a
  per-gene statistic; every gene is RANKED by how well its knockout recovers
  patients in the PCA/UMAP space built from the DF reactions (`mode='df'`).
  Relaxing the cutoff (26→116 reactions) REBUILDS that space: the PCA is refit on
  a different feature set → new axes, new geometry → "recovery toward control" is
  measured along different coordinates. RRA ranks are relative, so adding reactions
  can demote a former winner and promote a new one. SLC7A1 was the best lever in
  the transport-dominated 26-reaction space; once BCAA/FAO reactions enter at 0.3,
  DBT becomes the best lever and SLC7A1 falls away.

Analogy: relaxing the reaction cutoff is not "lowering the pass mark on the same
exam" (monotone) — it is "changing the exam questions" (the targets are scored in a
space the cutoff redefines).

Implication: the target list is a function of `--flux_rbc`. ALWAYS report the
cutoff; do NOT assume a 0.5 target reappears at 0.3 or vice versa. The rbc=0.3 CD
result is the more trustworthy of the two (116 reactions give a more stable
embedding, and the winners form a coherent BCAA/β-oxidation cluster with recurrence
0.61, bootstrap 1.0), but robustness still requires a cutoff-sensitivity analysis.
Two mitigations (see FUTURE_REQUESTS.md): (1) run the ranking at 0.3/0.4/0.5 and
keep only cutoff-stable targets; (2) optionally rank in the `mode='all'` reaction
space, which is cutoff-independent by construction.

---

## Session: 2026-06 (pipeline extension, DE-flux relaxation, CD vs UC target ranking)

### What was done
- Explored the full pipeline (`__main__.py` and modules) and documented it.
- Added DE-flux CLI flags (`--flux_pval`, `--flux_filter`, `--flux_rbc`) so the
  previously hard-coded thresholds (p<0.05, raw_p, rbc=0.5) are tunable.
- Added the cohort target-ranking module
  (`statistical_analysis/cohort_target_ranking.py`) + flags
  (`--target_fdr`, `--target_effect_min`, `--target_topk`). Outputs
  `cohort_target_ranking.csv` and `viz/moma_gene_scores_full.csv`.
- Fixed per-patient score normalization (raw MOMA scores were ~1e9–1e11 and
  saturated the statistics; now normalized inside the ranking step).
- Ran CD-ileum and UC-rectum at `--flux_rbc 0.3`.

### Runs / counts (current)
| Run | DE-flux reactions | n (test) | MOMA targets (significant) |
|---|---|---|---|
| CD-ileum, rbc=0.5 | 26 | 36 | — |
| CD-ileum, rbc=0.3 | 116 | 36 | **35** |
| UC-rectum, rbc=0.5 | 159 | 20 | — |
| UC-rectum, rbc=0.3 | 369 | 20 | **0** |
| (testable reactions: CD 2219, UC 2170) | | | |

### Key findings
- **CD produces a strong, coherent target signal; UC produces none.**
  CD top targets cluster on **mitochondrial BCAA catabolism (DBT, BCKDHA, BCKDHB)
  and fatty-acid β-oxidation (HADHA, HADHB, ECHS1, ACAA2, ACADM)** —
  recurrence up to 61% of patients, bootstrap 1.0, rra_fdr ≈ 1e-26 to 1e-54.
  Top-50 gene overlap between CD and UC = 0.

- **UC has MORE DE flux than CD at every cutoff (159 vs 26 at 0.5; 369 vs 116 at
  0.3), yet 0 targets.** Reasons established:
  - DE-flux count ≠ signal strength. UC's reactions are diffuse across many
    subsystems and patient-specific.
  - It is NOT a GPR artifact: UC has more gene-linked DE reactions than CD
    (289 vs 104) and more FAO DE reactions (21 vs 2).
  - UC's best gene recurs in only **2 of 20 patients** (recurrence_prop max 0.1),
    bootstrap ≤ 0.36 — genuinely non-reproducible, not just a borderline FDR.
  - Partly statistical power (n=20 UC vs n=36 CD) + low cohort agreement
    (Kendall W ≈ 0.17 both).

- **The CD/UC difference is NOT caused by the `target_switch_n` method switch.**
  `mixedlm_estimate` is 100% NaN in BOTH cohorts — the mixed-effects model never
  ran for either (CD at n=36 should have triggered it but did not / fell back to
  RRA). Both cohorts were scored by the SAME method (RRA + permutation).
  Applying CD's exact rule to UC still gives 0 (UC best rra_fdr = 0.18).
  → `--target_method auto` is currently inert; always lands on RRA.

- **Normalization fix verified working.** After the fix, `median_score` is on a
  sane scale (CD −7e-10..7.5; UC −11.9..1.2e-9) and `rank_biserial` is no longer
  saturated at ±1. The raw `moma_gene_scores_full.csv` still shows ~1e9 because
  normalization is applied inside the ranking computation, not written back.

- **`cohort_target_ranking` is derived from DE flux, not all flux** (MOMA recovery
  scored in the `mode='df'` PCA/UMAP space). Confirmed in `__main__.py` lines
  ~316–352. Therefore the ranking is cutoff-dependent.

### Known issues surfaced
> Update (2026-06-24): most of these are now resolved — see the 2026-06-24
> session above for status. (No separate FUTURE_REQUESTS.md exists; tracking
> lives here.)
- `--target_method` / `--target_switch_n` flags were NOT implemented (caused
  "unrecognized arguments" error on a UC run). — **resolved**, implemented
  with `--target_nperm` / `--target_nboot`.
- Mixed-effects model never fires (all-NaN) — convergence or guard issue.
  — **explained**: by-design guards (`switch_n`=40 > n; G > `mixedlm_max_genes`).
- RRA ranks genes with ~0 recurrence (RRA-vs-recurrence disconnect).
  — **mitigated** (recurrence reported + gated); raw RRA order can still surface
  low-recurrence genes.
- Option B (effect size on both test branches) not implemented. — **resolved**,
  Cohen's d added to the parametric branch with `cohens_d_cutoff`.

### Caveats
- UC analysis was run by the user; CD vs UC comparison treats both current
  rbc=0.3 runs. Earlier "UC=0 at 0.5" referred to the OLD bundled example data,
  not the current IBDMDB rectum inputs — do not mix those numbers.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               