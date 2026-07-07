# FUTURE_REQUESTS — tf_metabolism

Known issues, planned improvements, and design decisions still open.
Priority: [P0] blocker / [P1] high / [P2] medium / [P3] nice-to-have.

---

## [P1] Multi-model / multi-MEM — Phase-2 follow-ups
Phase 1 (done) added `--model {Recon2M.2,Human1,Human2}` and `--mem_method`
(multi-value), a per-model resource registry (`metabolic_model/model_config.py`),
Human-GEM cloned at pinned tags (Human1=v1.0.0/e4feea4, Human2=v2.0.0/635f533),
model-inspection reports, Entrez↔Ensembl omics mapping from the in-repo
`genes.tsv` (Human2 82.8% / Human1 83.5% model-gene coverage), output tagging
`<-o>_<model>_<mem>`, and fail-loud on unsupported MEMs. Recon2M.2 + tINIT is the
default and unchanged. Remaining:

- **[P1] Human-GEM tINIT reconstruction — RAVEN checkTasks port DONE; blocked now on a growth medium.**
  Phase-1b/c (DONE): four resource adapters, namespace resolver, validate-then-run gate.
  Phase-1d (DONE): a faithful **RAVEN `checkTasks` engine** (`evaluate_human_tasks` in
  `model_config.py`) that handles the constructs the old internal format could not —
  **EQU** (inject the forced reaction per task with its bounds), **simultaneous
  multi-OUT** (force all output demands together in one feasibility LP), **ALLMETSIN**
  (open all uptake/secretion), **SHOULD-FAIL** (pass iff infeasible), plus a base
  inorganic + amino-acid medium (matching the Recon `MetabolicTasks_bigg.csv` Medium).
  `validate_raven_conversion` runs it as the cheap oracle on the FULL model.
  **Verified: all 57 Human2 essential tasks resolve and PASS on the complete model,
  0 fail, in ~40 s** (`RAVEN_CONVERSION_VALIDATION.txt`, overall=PASS). Fail-loud
  preserved: any unresolved metabolite/EQU token → the evaluator raises → gate FAILS.

  Phase-1e (DONE): growth-medium / biomass-feasibility adapter. **Diagnosis:** the
  Human-GEM biomass reaction `MAR04413` loads **bounded to [0, 0]** (standard GEM
  convention of disabling the objective), while all 1660 exchanges load **open for
  uptake** (Human-GEM's generic default medium) and all 38 biomass precursors are
  producible — so biomass was 0 purely because the objective reaction was shut, not a
  medium gap. **Fix (`apply_human_medium`):** open `MAR04413`; keep Human-GEM's shipped
  exchange bounds (no ad-hoc medium, no gap-fill). `biomass_feasibility` is wired into
  the pre-flight as a blocking growth oracle. **Verified end-to-end: full Human2 model
  biomass = 314; the gate un-gates and a capped tINIT smoke reconstructs one sample
  (C3002) → a real context model, 9808/12931 reactions (subset), 2848 genes, biomass =
  241 (growth-feasible), r.subsystem populated on all 9808** → `functional_C3002.xml`.
  **Human2 + tINIT now reconstructs.** Remaining:
    1. Full-cohort reconstruction + downstream — **DONE (Phase 1).** `--model Human2
       --mem_method tINIT` now runs the SAME end-to-end pipeline as Recon2M.2 into the
       tagged folder: full per-sample reconstruction (`reconstruct_human`, all samples,
       skip-guarded, per-model sanity gate biomass>0 & objective==MAR04413) → flux
       prediction → DE flux → flux-sum → gene-subsystem-flux → pathway enrichment (over the
       138 attached Human-GEM subsystems) → MOMA (clean ~2,184 ENSG gene set) → cohort
       ranking → viz. Three Human-only adaptations, all guarded so Recon2M.2 + tINIT stay
       byte-identical: (a) `flux_prediction.calculate_flux`/`predict_metabolic_fluxes`
       gained a `biomass_reaction_id` param (default `'biomass_reaction'`) + tolerate
       `getpra_file=None` (GeTPRA is Recon/BIGG-only); (b) subsystems attached onto the base
       model for enrichment/flux-sum; (c) reconstruction routed to `reconstruct_human`.
       A fail-loud namespace/objective oracle at each stage boundary writes
       `HUMAN2_DOWNSTREAM_VERIFY.txt` (reactions MAR####, genes ENSG, MOMA objective
       MAR04413). Verified end-to-end on a 2v2 CD-ileum subset (reconstruction/flux/flux-sum/
       enrichment/MOMA all run; DE flux needs real cohort n for significance). Full cohort is
       Gurobi-heavy (MOMA ≈ n_condition2 × ~2,184 gene-KOs).
    2. Human1 v1.0.0 ships no essential-task `.txt` (only `LeakCheck.xls`) — extract its
       tasks from the `.mat`/`.xls` for parity, then the Human1-vs-Human2 comparison run.
  Also: Human1 v1.0.0 ships no essential-task `.txt` (only `LeakCheck.xls`) — extract
  its tasks from the `.mat`/`.xls` to bring it to parity, then the Human1-vs-Human2 run.

- **Model-counting parity (diagnosed; see MODEL_COUNTING_PARITY.txt).**
    - *Human context orphan metabolites AND genes — FIXED.* The Human smoke path left
      reaction-pruned but metabolite-/gene-unpruned models. cobra 0.5.4 keeps BOTH orphan
      metabolites and orphan genes in memory after reactions are removed, so the context
      reported the BASE counts (8,461 mets / 2,848 genes). Fixed in `reconstruct_human_smoke`
      (Human-only) by reusing the Recon convention exactly: the manual orphan-metabolite loop
      (`tINIT.reconstruct_GEMs` :131-133) PLUS a write→reload (:154-155) so cobra rebuilds
      `model.genes` from surviving GPRs (with `use_fbc_package=False`, orphan genes are never
      written). Verified across two cohorts, 0 orphans of either kind: **CD-ileum C3022 =
      9,563 rxns / 2,184 genes / 6,700 mets**; **UC-rectum C3022 = 9,253 / 2,157 / 6,598**
      (genes drop from 2,848 → true active ~2,160–2,185; mets 8,461 → ~6,600–6,700). Recon +
      tINIT byte-identical.
    - *Round-trip objective regression — FIXED.* `write_sbml_model(use_fbc_package=False)`
      round-tripped Human-GEM's SHIPPED default objective (`MAR13082`) into the file, so the
      on-disk model downstream reads optimized the wrong reaction. Fixed by setting
      `ctx.objective = MAR04413` (and keeping it open) BEFORE the write. Verified on the
      fresh-read file: objective=`MAR04413`, bounds `[0,1000]`, biomass FBA=242, and the
      138 subsystems survive the round-trip (100% reaction coverage). Because
      `reconstruct_human_smoke` returns FILE PATHS (not the in-memory object), any downstream
      (MOMA, future WES/methylation constraints) reads this clean reloaded model — clean gene
      set, re-attached subsystems, correct biomass objective. NOTE: no Human downstream is
      wired yet, and there is currently no WES/methylation constraint code in the repo
      (transcriptomics-only) — the guarantee here is structural (disk read), to be re-checked
      when that downstream is built.
    - *Recon2M.2 load is BIGG-core — INTENTIONAL, documented.* The pipeline loads the
      BIGG-namespace Recon2M.2 (`Recon2M.2_Entrez_Gene_BIGG.xml`: 5,842 rxns / 3,368 mets
      / 1,663 genes). All five shipped `*_BIGG`/`*Transcript` variants share the same
      5,842/3,368 (only gene namespace differs), so it is not a truncated load. The ~25%
      reaction / ~37% metabolite gap vs canonical Recon2M.2 (7,785 / 5,324; Ryu 2017) is
      BIGG-mapping loss — features without BIGG IDs dropped in conversion — and is required
      because the pipeline's tasks/medium/metabolite lists are all BIGG-named. Not a bug;
      do not change the Recon load.
- **[P2] Real multi-MEM backend (PART B).** Only tINIT runs here (cobra 0.5.4 /
  py3.6). Build a SEPARATE modern-Python env with `troppo` (GIMME/iMAT/FASTCORE/
  (t)INIT) bridged via subprocess; MBA/mCADRE additionally need MATLAB + COBRA
  Toolbox. Until then non-tINIT `--mem_method` values fail loud by design.

---

## [P0] MOMA target score must measure RECOVERY, not disruption magnitude
**Problem (validity blocker).** The cohort target ranking scores genes by
`recovery_proj_pca` on a raw ~1e11 scale, which is dominated by the MAGNITUDE of
the flux perturbation, not by whether the knockout moves CD patients toward the
nonIBD control state. Consequence: the CD "top targets" (DBT, BCKDHA, BCKDHB,
HADHA/B, ECHS1, ACADM) are essential/high-flux genes whose KO causes a big
disruption — by the honest distance metric they barely recover patients (DBT:
distance-to-control 2.85→2.70, only 53% of patients closer; KO-shift UMAP shows
the cohort collapsing AWAY from control). **These targets are NOT validated; do
not carry them forward.** See WORKLOG 2026-06-24b.

**Literature note (which metric is conventional).** MOMA itself (Segrè 2002) only
defines the WT↔KO flux distance — it has NO "toward control" notion; that is an
add-on analysis layer. In the KO-simulation reversal / signature-reversal field,
the conventional recovery metric is a **normalized DIRECTIONAL projection onto the
disease→control axis** (Signature Projection Score; TReD), which is explicitly
preferred over raw Euclidean distance because it captures direction AND magnitude
relative to control. So `recovery_proj_pca`'s *concept* (a directional projection)
is the standard — the bug is only that this pipeline's version is UNNORMALIZED
(~1e11), collapsing it into a magnitude score that even disagrees with actual
proximity. Cancer GEM KO papers (PNAS 2025 breast-cancer KO+UMAP+k-means; redox
radiation-resistant tumors) also judge recovery by cluster/region membership, not
raw distance. Refs in WORKLOG 2026-06-24c.

**FINAL FIX — competitive null (MTA-style score + Reporter-Metabolite-style
background restandardization). This SUPERSEDES the intermediate fixes below.**
The decisive problem is a GLOBAL DRIFT: almost any knockout nudges most patients
slightly toward control (cutting flux pulls the hyperactive disease state toward
the lower-flux control region). So ALL absolute tests are confounded — "recovery
> 0" flags ~everything; the `frac_closer > 0.5` / `proj_vs_distance_agree` majority
gate is knife-edge at n=36 (DBT 0.472 vs SDSL 0.528 = 2 patients, both at coin-flip,
binomial P≈0.4–0.7); and ranking by mean recovery is outlier-driven. Verified: the
"significant" genes are NOT above a coin flip.
The correct, literature-standard test is a **competitive null**:
1. **Score** each knockout for movement toward the control centroid (MTA/rMTA-style,
   directional, in 30-PC PCA space).
2. **Restandardize against a random/sampled-knockout background** (Reporter-
   Metabolite-style: subtract mean, divide by SD of the background of random
   knockouts of matched perturbation magnitude), then FDR-correct. A gene is a
   target only if it recovers MORE than a random knockout — which removes the drift.
Refs: Patil & Nielsen 2005 PNAS (reporter metabolites / random-set background);
Yizhak MTA 2013 & Valcárcel rMTA 2019 (disease→healthy KO ranking, relative);
competitive-vs-self-contained null (gene-set theory). See WORKLOG 2026-06-24d.
**Prototype result (drift-controlled permutation null on existing data): 0/1362 CD
genes pass p<0.05 or FDR<0.05; DBT/FASN/SDSL/BCAA-FAO all p≈0.77.** Honest, final
CD conclusion: **no single-gene knockout robustly reverts the CD metabolic state
beyond a non-specific perturbation effect.** Implement the competitive null as the
official `significant` rule; expect empty CD target set.
- No MOMA/Gurobi rerun needed; regenerate from existing `UMAP_Results_*.csv`.

**Superseded intermediate attempts (kept for history — do NOT reuse as the test):**
- raw `recovery_proj_pca` magnitude score (rewarded disruption, not recovery);
- normalized directional projection alone (better, but still absolute → drift-
  confounded);
- `frac_closer>0.5` / `proj_vs_distance_agree` majority gate (knife-edge, unsound
  at n=36);
- mean-recovery ranking (single-patient-outlier driven).
All produced apparent targets (BCAA/FAO; then DHFR/FASN/folate) that DO NOT survive
the competitive null. These shortlists are artifacts — do not carry forward.

## [P0b] KO plots collapse — fit_transform on baselines+KOs (display only, not significance)
**Problem.** Both KO-plotting functions place knockout points with
`reducer.transform()` into a baseline-ONLY UMAP (`run_moma_results_umap` L536, per-
patient; `plot_gene_ko_shift` L1254, ko_shift). Knockout fluxes are out-of-
distribution, so transform() collapses them to a corner (active genes) / WT pile
(silent genes, ~98%) — gene-independent (FASN==DBT). It is a visualization artifact,
NOT UMAP-vs-PCA and NOT biology. The baseline UMAP (`run_r_umap` L328) uses
`fit_transform` and is fine. The reference IBD project's FASN scatter works because
it fit the UMAP on the KO cloud (fit_transform).
**Fix.** In BOTH KO-plotting functions, build the embedding with **fit_transform on
baselines + KO points together** (KOs in-distribution → scatter), instead of
transform() into the baseline-only model. Accept that baseline dots will not match
the existing `UMAP_df` coordinates (standalone KO plot). Optionally also offer a
recovery-axis plot (x = signed recovery toward control, y = orthogonal residual, in
PCA) where horizontal length literally = how much closer to control.
**Boundary (critical).** This is a DISPLAY fix only. Significance is decided by the
competitive null ([P0]) on the high-D recovery score — NEVER by 2D proximity in any
embedding. Keep display and significance strictly separate.
**Scope.** No MOMA/Gurobi rerun; regenerate plots from existing
`targeting_results/` + flux. Remove superseded absolute-recovery gating
(frac_closer/proj_vs_distance_agree/mean-ranking) — see [P0].

**Contingency — PRESERVE the removed stat analysis (do not lose it).** The
absolute-recovery code being removed (`frac_closer`, `frac_delta_pos`,
`proj_vs_distance_agree`, mean-recovery ranking, and the normalized directional
projection score) must be kept recoverable in case the `fit_transform` display fix
does NOT resolve the collapse and we need to fall back. Keep it by EITHER (a) moving
it into a clearly-named deprecated module/function (e.g.
`cohort_target_ranking_legacy.py` or a `# DEPRECATED — fallback only` block) rather
than deleting outright, OR (b) ensuring it lives in a committed git revision/branch
referenced here. Reference implementation snippet (the gate logic) for restoration:
`frac_closer = mean(recovery_proj_norm > 0)`; `frac_delta_pos = mean(delta_distance
> 0)`; `proj_vs_distance_agree = (frac_closer>0.5) & (frac_delta_pos>0.5)`;
`significant_legacy = rra_fdr<fdr_cutoff & recurrence_prop>=effect_min & frac_closer>0.5
& proj_vs_distance_agree`. These remain SUPERSEDED for significance (see [P0]) — kept
only as a fallback if the display fix fails, never as the primary target test.

## [P1] Add `--target_method` / `--target_switch_n` flags — DONE (verify)
Implemented in the 2026-06-24 session (auto/rra/mixedlm + switch_n/nperm/nboot).
Leave here until confirmed on a fresh run.

## [P1] Fix the mixed-effects model never running
`mixedlm_estimate` is all-NaN on CD/UC: (a) switch_n=40 > CD n=36, and (b) G≈1362
genes > mixedlm_max_genes=300, so it is skipped. By design for now, but if a
mixed-effects view is wanted at large n, make the gene set tractable (e.g. restrict
to candidate genes) or lower switch_n with care.

## [P1] Cutoff-sensitivity analysis for the target ranking
The ranking is built in the DE-flux (`mode='df'`) space, so it is cutoff-dependent
and NON-monotone: SLC7A1 was significant at rbc=0.5 but not 0.3; DBT the reverse
(WORKLOG note "Target ranking is cutoff-dependent and NON-monotone"). Run the
ranking at rbc 0.3/0.4/0.5 and keep only cutoff-stable targets. Alternatively add a
`mode='all'` ranking that is cutoff-independent by construction.

## [P1] Decouple the effect-size gate from the test choice (Option B) — DONE (verify)
Implemented: Cohen's d on the parametric branch with `cohens_d_cutoff`, alongside
rank-biserial / `rbc_cutoff`. DE-flux call passes `cohens_d_cutoff=0.8`. Verify on a
large-n cohort that the effect-size gate stays active at n≥30.

## [P2] Resolve the RRA-vs-recurrence disconnect
RRA can rank genes that are top-k in almost no patient. Mitigated by reporting and
gating on recurrence; consider requiring a minimum effect magnitude before a gene
enters the RRA list.

## [P2] Parallelize the per-patient MOMA loop
`run_targeting_simulation` is O(patients × genes × LP) and serial — fine for n≤36,
slow for large cohorts. Embarrassingly parallel (joblib / Snakemake). Statistics
unaffected.

## [P3] Tagged / non-overwrite outputs for multiple cutoffs
Output filenames are fixed, so a new cutoff overwrites prior results. Workaround =
separate `output_*_rbc_*` folders pre-seeded with the expensive inputs. Optional
`--tag` flag would let cutoffs coexist (needs threading through the R figure scripts).

---

## Open scientific questions

- **Are there ANY robust CD targets?** Prototype says no gene recovers >72% of
  patients; re-evaluate after the [P0] score fix. Validate DHFR/FASN/folate against
  IBD literature if they survive.
- **UC null:** real diffuse biology vs underpowered cohort (n=20, Kendall W≈0.17)?
- **CD vs UC contrast** must be re-examined after [P0] — the previous "CD strong,
  UC null" story was partly a magnitude artifact.

---

## Standing reminders
- DE flux is raw_p + effect size, NOT FDR — all DE-flux and targets are exploratory.
- `cohort_target_ranking` depends on the DE-flux set and the `--flux_rbc` cutoff.
- Target score (post-[P0]) = signed recovery toward control, not disruption magnitude.
- Outputs are gitignored and not recoverable from git — back up before deleting.
- Do not mix old bundled-example numbers with current IBDMDB numbers.
