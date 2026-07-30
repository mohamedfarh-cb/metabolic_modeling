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
### STATUS 2026-07-29: PARTIALLY RESOLVED — stages 1-2 done, the FINAL FIX (competitive null) is NOT implemented
Verified by re-reading the code and recomputing from the existing Human2 runs
(see WORKLOG 2026-07-29). **Done:**
- `recovery_proj_pca` is now a signed projection onto the **unit** Test→Control
  centroid axis in PCA space (`r_based_viz/moma_project.R:119-120`, applied :208).
  Scale is sane: CD 4.81..23.8, UC −32.2..20.5 — not ~1e11.
- The condemned absolute gates (`frac_closer>0.5`, `proj_vs_distance_agree`,
  mean-recovery ranking) are GONE. `significant` = `rra_fdr < cutoff AND
  recurrence_prop >= effect_min` only.
- Two partial drift controls exist: `_robust_scale_per_patient` centers each
  patient on its MEDIAN knockout ("recovers more than the typical KO"), and RRA
  runs on within-patient RANKS — inherently competitive against that patient's
  all-gene background, and invariant to the common positive offset.
- The old failing shortlist (DBT/BCKDHA/BCKDHB/HADHA/B/ECHS1/ACADM) no longer
  appears; Human2 CD targets are glycolysis/pyrimidine/transport, UC an oxidase +
  FAO arm, with 85-100% of patients moved toward control.

**NOT done — the FINAL FIX below still stands.** The background is ALL genes, not
a **magnitude-matched** random-KO set, so "recovers a lot because it is a large
perturbation" is still not separated from "recovers more than an equally large
random perturbation". Fresh evidence (2026-07-29, patient C3002, n=2173 genes):
global drift is real (ALL 2285 genes × 36 CD patients project positive, min
+4.81); cohort-wide Spearman(‖KO−WT‖, recovery) is only **0.112** — so recovery is
NOT merely magnitude — **but every one of the top-10 recovery genes sits at the
95th-98th percentile of perturbation magnitude** (PKM/PKLR 97.3, ENO1/2/3 98.1).
Also still open: the prototype's honest conclusion (0/1362 CD genes pass a
drift-controlled null) has NEVER been re-tested against the CURRENT scoring.
**Until the competitive null runs, the Human2 target lists remain
hypothesis-only — do not report them as validated targets.**

**EVIDENCE 2026-07-29c/30 — the drift is real, confirmed on a published metric.**
`prototype_rmta_ts.py` implements rMTA's Transformation Score (Valcárcel 2019 Eq. 8;
bounded [−1,+1], normalized by total change, so magnitude cannot win) on the existing MOMA
output. Result: **CD background median mTS +0.2563 with 2284/2285 genes (100%) positive;
UC −0.0327 with 55/2258 (2%) positive.** So the drift is a property of the CD run, not of
any gene, and CD's raw scores (~18× UC's) are uninterpretable without their background —
UC's excess over background is ~3× CD's. Independent corroboration: the top of the filtered
UC list is UC's own oncoprint oxidase arm (ACOX1 +0.090, ACOX3 +0.089, NOS2 +0.087,
DAO/PAOX +0.085), recovered by a metric that never saw the RRA ranking.

**Terminology correction:** "global drift" is THIS REPO's wording, not a literature term.
The published concern is rMTA's *under-determination* (top MTA knockouts move as readily in
the undesired direction), and the field's fix is a **worst-case control (wTS)** — not the
magnitude-matched random-KO background this item prescribes (that is the CMap-tau
approach). Both are defensible; cite whichever you implement. The prototype covers only
rMTA's **mTS**; bTS/wTS need new MTA MILP runs per gene per patient and are NOT derivable
from stored MOMA output. **[P0] therefore remains OPEN.**

**Original problem (historical — the ~1e11 score is fixed; kept for context).**
The cohort target ranking scored genes by
`recovery_proj_pca` on a raw ~1e11 scale, which was dominated by the MAGNITUDE of
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
  confounded) — **NOTE: this is what the pipeline CURRENTLY implements** (plus
  per-patient median-centering and rank-based RRA, which remove the drift OFFSET
  but not the magnitude confound). It is stage 2 of 3, not the final test;
- `frac_closer>0.5` / `proj_vs_distance_agree` majority gate (knife-edge, unsound
  at n=36);
- mean-recovery ranking (single-patient-outlier driven).
All produced apparent targets (BCAA/FAO; then DHFR/FASN/folate) that DO NOT survive
the competitive null. These shortlists are artifacts — do not carry forward.

## [P0c] The desired direction is a COHORT average applied to every patient — ~35% of each patient's reactions are on the wrong side of it
### FOUND 2026-07-31 by `prototype_control_distance.py`. Affects BOTH competing scores; not fixed.

`Differential_fluxes.csv Status` is a population-level statement ("up in disease **on
average**"). Every downstream score applies that single per-reaction direction to every
patient: the pipeline's `recovery_proj_pca`, and the rMTA mTS prototype alike. Measured
across all DE reactions x all patients:

| | patient on the DISEASE side (cohort direction applies) | patient already on the CONTROL side (cohort direction is WRONG for them) |
|---|---|---|
| CD-ileum (244 rxn x 36 pt) | 61.8% | **37.7%** |
| UC-rectum (771 rxn x 20 pt) | 62.0% | **33.8%** |

Per-patient the disease-side share ranges 47.5-75.8% (CD) and 45.4-71.1% (UC) — for the
worst CD patient the cohort direction is wrong for **52.5%** of their reactions.

**Consequence, measured.** Absolute distance to control (direction-free,
`sum |v_i - c_i| / gap_i`) **increases for every knockout in every patient** — CD median
+174.9%, UC +74.2% for a typical KO. Decomposing the 344 reactions that end further away
for a median-perturbation gene (UC C3003): **58% because the patient was already on the
control side**, 42% wrong-way movement or overshoot (59 reactions crossed past control).

So mTS can be positive for 100% of CD genes while the patient's actual distance to control
grows. Those are not contradictory — one measures agreement with a cohort-average
direction, the other measures where the patient ended up. **Neither the pipeline score nor
mTS can detect this, because both assume the cohort direction.**

**What to do (unresolved).** Options: (a) per-patient desired directions (each patient's
own deviation from the control mean) — patient-specific but noisier and no longer tied to
the DE test; (b) keep cohort directions but weight each reaction by how consistently
patients sit on the disease side; (c) report absolute distance alongside any directional
score so the disagreement is visible rather than assumed away. Do NOT quietly pick one —
each changes what "toward control" means.

Not a defect of the distance metric: verified no units mismatch (median |WT| 38.2 vs |KO|
26.9; perturbation 0.5x the WT-to-control gap; 0/2148 genes leave flux unchanged, which is
expected — ROOM notes MOMA redistributes ALL fluxes on any knockout).

## [P0b] KO plots collapse — fit_transform on baselines+KOs (display only, not significance)
### STATUS 2026-07-29: RESOLVED on the R path (the one all four runs used); Python fallback still unfixed — see the new [P2] below
Commit `e388856` ("MOMA plots: project reference set through predict() for one
consistent frame") put KO points, the WT cross and the reference baselines in ONE
frame: `moma_project.R:198-203` projects the WT through the SAME `predict()`
pipeline as the KOs (delta = 0), and `wt_x`/`wt_y` are written per sample so the
cross is drawn at the projected position rather than the training-fit position.
The displacement vector is therefore self-consistent. All four IBDMDB runs used
the R backend (their `UMAP_Results_*.csv` carry `recovery_proj_pca`), so the
delivered figures are fine. Original problem below, for context.

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

## [P1] Adopt the calibrated KO figure into the pipeline (prototype exists and works)
`prototype_ko_cohort_calibrated.R` + `prototype_rmta_ts.py` (both repo root, standalone,
uncommitted) demonstrate the fix for the §3i figure's biggest interpretation hazard: the
"k/n moved toward control" count reads as evidence when it is baseline behaviour. Verified
on all 9 ELOVL-set genes × 2 cohorts (WORKLOG 2026-07-30).

**Adopt, with these three changes — do not port the prototype verbatim:**
1. **Keep** the cohort-background clause in the subtitle (`[cohort background: 79% of ALL
   genes move toward control]`). Single highest-value addition.
2. **Keep** the inset histogram of the all-gene mTS distribution with the gene marked. It
   makes the competitive-null argument in one glance and needs no prose.
3. **DROP the percentile** (`beats 97% of genes`). It does NOT separate the cohorts —
   CD 95-97% vs UC 98-99% — so a reader quoting it concludes CD is nearly as good. The
   separation lives in the excess over background, not the rank.
4. **Cut the subtitle to 2 lines.** The prototype's is 5 and dominates the panel.

**Blocker on full adoption:** the mTS values come from `prototype_rmta_ts.py`, which reads
the full stored MOMA matrices (~5.8 GB across both cohorts, ~10 min per cohort). Wiring the
figure to require that file means either (a) running the mTS step inside `run_pipeline`
after MOMA, or (b) keeping the figure's mTS panel optional and rendering it only when
`viz/prototype_rmta_mTS.csv` exists. (b) is cheaper and keeps the flag opt-in.

Also fold in: a **min-patient filter is REQUIRED** in the mTS table (genes present in 1-2
patients' models otherwise top the ranking — CD DCT/GNMT/SARDH, UC ADSS2/TNRC6B/SLC17A7).
All reported numbers use >=50% of patients.

## [P2] Python MOMA-projection fallback never got the [P0]/[P0b] fixes — silent regression risk
`utils.run_moma_results_umap` (the non-`use_r` branch, `utils.py:642-657`) writes
only `X`, `Y`, `distance_to_center` — **no `recovery_proj_pca`, no WT reference
projected through the same frame**, and it still uses `reducer.transform()` into
the baseline-only model (the [P0b] artifact). `export_moma_gene_scores_full`
(`utils.py:885-892`) then falls back to `_l2_score_from_moma` — the raw
displacement magnitude, i.e. exactly the disruption-not-recovery score [P0]
exists to eliminate — and does so with only a `logging.warning`, so a run would
silently produce a ranking built on the discredited metric.

No delivered result is affected: all four IBDMDB runs used the R backend
(verified — their `UMAP_Results_*.csv` carry `recovery_proj_pca`). This is a
latent trap for any future R-less run.

**Fix options:** (a) port the recovery-axis computation (unit Test→Control PCA
vector + WT-through-the-same-transform) into the Python path so both backends
agree; or (b) if the R path is the only supported one, make the fallback FAIL
LOUD instead of silently degrading to the L2 score. (b) is cheap and matches the
repo's fail-loud convention elsewhere; (a) is correct if the Python path is meant
to be a real backend.

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

## [P2] Cross-model compartment comparison must normalize to background (not raw %)
The Recon2M.2-vs-Human2 compartment figure
(`analysis_compartment_recon_vs_human2/`) currently compares the RAW organelle %
of DE-flux genes / targets. Raw % is NOT comparable across models, because each
model has a different whole-model compartment background and a different gene
inventory:
- Whole-model mitochondrial background: Recon2M.2 ≈19% vs Human2 ≈12%. Human2 is
  ~2× larger, but the genes it adds over Recon are disproportionately
  cytosolic/nuclear/other, so every downstream mito fraction starts from a lower
  floor. This alone flips the naive reading: UC DE-flux is 35% mito in Recon vs
  15% in Human2, yet in ABSOLUTE counts Recon has MORE mito DE-flux genes
  (143 vs 81) — Human2's extra DE-flux genes went to cytosol.
- Fix: report ENRICHMENT = DE-flux organelle% ÷ whole-model organelle% (fold over
  each model's own background). Normalized, both models AGREE UC DE-flux is
  mito-enriched (Recon ×1.78, Human2 ×1.24) and cytosol-led; the raw bars oversell
  the divergence. Regenerate the figure as background-normalized fold-change bars
  for any cross-model claim.
- Deeper caveat: the histograms summarize INDEPENDENT DE-flux reaction sets
  (BIGG vs MAR — the models don't share reaction ids), and each gene is collapsed
  to a single dominant compartment (mode over its reactions). Prefer a
  reaction-level compartment comparison, and map to a shared id space (MetanetX)
  before comparing which reactions/pathways — not just which compartment % — agree.
- Reliable across models = DIRECTION (cytosol-led both diseases; UC mito-enriched;
  CD not). NOT reliable = the exact percentages.

### Pathway-name harmonization (blocks the enriched-pathway comparison)
The comparative enriched-pathway barplots
(`analysis_compartment_recon_vs_human2/enriched_pathways_barplot_{CD,UC}.png`)
match pathways by EXACT subsystem name, but the two models use different
subsystem vocabularies, so the cross-model overlap is deflated by synonyms, not
just biology. At FDR<0.1 only 1 pathway (CD) / 0 (UC) matched in both — largely a
naming artifact. Known synonym pairs:
- Recon `Transport, extracellular` ≈ Human2 `Transport reactions`
- Recon `Beta oxidation of even-chain fatty acids` ≈ Human2 `Fatty acid oxidation`
- Recon `Valine, leucine, and isoleucine metabolism` → no exact Human2 twin
Fix options: (a) a manual Recon↔Human2 subsystem synonym map to mark truly shared
pathways in the barplot; (b) better, map reactions to a shared ontology (MetanetX /
KEGG / Reactome) and re-derive enrichment on the shared pathway space so overlap is
computed on identity, not strings. Until then, read the barplots for magnitude/
direction WITHIN each model, and treat the low cross-model overlap COUNT as
inflated-negative. The genuinely Human2-only lipid biology (bile acid, omega-3 /
unsaturated FA, peroxisomal carnitine shuttle in UC) is real — Recon's vocabulary
cannot express it — and survives this caveat.

### Concordance analyses added (context for the above)
`analysis_compartment_recon_vs_human2/` now holds cross-model comparisons:
`deflux_gene_overlap` (Jaccard CD 0.13 / UC 0.18), `deflux_direction_concordance`
(shared-gene sign agreement CD 67% ρ=0.39 > UC 58% ρ=0.20 — CD signature more
model-robust), and `moma_target_rank_concordance` (Spearman ρ≈0, 0/20 consensus —
MOMA target RANKS do NOT reproduce across GEMs; single-model target lists are
hypothesis-only). Plus `model_size_boxplot_by_subtype` (Human2 ~2.8× reactions,
1.7× genes, 3.1× metabolites). Folder is NOT gitignored.

---

## [P2] metabolic_related_rna_seq.R flux violins — selection + interpretation caveats
The per-run flux-violin scripts (CLAUDE §3h) accumulated several deliberate design choices
worth revisiting:
- **"Owner-by-region" gene selection is a heuristic.** Each pathway is sourced from the
  cohort(s) with ≥1 DE-flux oncoprint target in that subsystem (glycolysis→CD, omega-3→UC,
  sphingolipid→both), then the same set is shown in both. A pathway with no DE-flux target in
  either cohort falls back to "both". This encodes the CD/UC region split but is not derived
  from a formal region classifier — if the oncoprint targets shift (e.g. after the [P0] score
  fix) the ownership can flip. Consider deriving ownership from `uc_cd_metabolic_comparision`
  region output instead of recomputing inline.
- **Sphingolipid is NOT significantly enriched in UC** (Up FDR 0.231; CD 0.098) — it is shown
  on the UC side only because SLC33A1 is a shared DE-flux target. Decide whether to keep it on
  the UC panel or mark it "target-driven, not pathway-significant".
- **Between-violin bracket ≠ the `*`.** `*` = subtype-specific oncoprint target; the bracket =
  per-gene Wilcoxon on the plotted (subsystem-averaged) flux. They diverge because `de_flux`/
  target status is reaction-level across ALL subsystems while the plot averages ONE subsystem's
  reactions (PKM, HADHA, ACOX3 are starred with no bracket). This is intended but easy to
  misread — the subtitle must keep both definitions explicit.
- **Not-plottable target genes** are dropped + logged (e.g. omega-3 ELOVL3/HACD1-4/TECRL,
  glycolysis PGAM2/PGAM4) — their reactions sit in a different subsystem than the ranking's
  `subsystem` label, so there is no in-subsystem flux to plot.
- **FASN / de-novo FA synthesis is absent** from the omega-3 panel (it is filed under
  `Fatty acid biosynthesis`). If de-novo synthesis matters, add a `Fatty acid biosynthesis`
  violin (would include FASN — DE-flux in UC — SCD, ACACA, etc.).

---

## Open scientific questions

- **Are there ANY robust CD targets?** Prototype says no gene recovers >72% of
  patients; re-evaluate after the [P0] score fix. Validate DHFR/FASN/folate against
  IBD literature if they survive. *(2026-07-29: still open. The scoring HAS changed
  — CD now yields 21 significant genes, glycolysis/pyrimidine/transport — but until
  the competitive null runs these are not the answer to this question.)*
- **UC null:** real diffuse biology vs underpowered cohort (n=20, Kendall W≈0.17)?
  *(2026-07-29: note the old "UC = 0 targets" statement is stale — the Human2 UC run
  yields 40 significant genes under the current scoring. The Recon-era UC null and
  the Human2 UC result must not be quoted interchangeably.)*
- **CD vs UC contrast** must be re-examined after [P0] — the previous "CD strong,
  UC null" story was partly a magnitude artifact. *(2026-07-29: partly re-examined —
  under Human2 + the current score both cohorts yield targets with near-zero overlap
  (shared = SLC33A1 only), so the original "CD strong / UC null" framing no longer
  holds. Final re-examination still waits on the competitive null.)*

---

## Standing reminders
- DE flux is raw_p + effect size, NOT FDR — all DE-flux and targets are exploratory.
- `cohort_target_ranking` depends on the DE-flux set and the `--flux_rbc` cutoff.
- Target score (post-[P0]) = signed recovery toward control, not disruption magnitude.
- Outputs are gitignored and not recoverable from git — back up before deleting.
- Do not mix old bundled-example numbers with current IBDMDB numbers.
