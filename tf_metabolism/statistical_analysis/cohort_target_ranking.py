"""Cross-cohort MOMA target ranking.

Size-adaptive design & the effect-size gate
--------------------------------------------
Per-patient MOMA single-gene-knockout simulations each yield a ranking of
candidate metabolic target genes (scored by how strongly a knockout pushes the
patient flux state toward the Control center).  This module aggregates those
per-patient rankings into ONE cohort-level ranking whose statistics stay valid
from small (n < 30) to very large (n >> 100) cohorts.

The guiding principle mirrors
``statistical_comparison.two_grouped_data_comparison``, which auto-switches
parametric/non-parametric inference at n = 30:

  * n_patients < 30  -> rank-based inference is primary: Robust Rank Aggregation
    (RRA) plus a within-patient label-permutation null.  Both are distribution
    free and robust when the cohort is underpowered.
  * n_patients >= 30 -> ADDITIONALLY fit a linear mixed-effects model (per-gene
    fixed effects, random patient intercept) for per-gene estimates.  RRA and the
    permutation null are still reported at every size, and the mixed model falls
    back gracefully to RRA-only if it fails to converge.

The effect-size gate is what keeps results meaningful at large n, where p-values
collapse toward 0 for trivially small effects.  A gene is flagged ``significant``
ONLY when it passes BOTH an FDR threshold AND a minimum effect size (the
recurrence proportion).  No statistic is ever reported as significant without an
accompanying effect size.

GPR non-independence
--------------------
Genes that share a reaction through the gene-protein-reaction rule are perturbed
identically by MOMA (e.g. PYCR1/PYCR2) and therefore receive identical scores in
every patient.  Such genes are collinear, not independent evidence; we detect
identical-score gene groups and tag them with a shared ``gpr_group`` so they are
reported together rather than double-counted.  Their individual p-values are NOT
independent and should be read at the group level.

Reference
---------
Kolde R, Laur S, Adler P, Vilo J. "Robust rank aggregation for gene list
integration and meta-analysis." Bioinformatics. 2012;28(4):573-580.
"""

import os
import logging

import numpy as np
import pandas as pd
from scipy.stats import beta, binom, wilcoxon, rankdata
from statsmodels.stats.multitest import fdrcorrection


# ---------------------------------------------------------------------------
# Cohort-level statistics (helpers)
# ---------------------------------------------------------------------------
def _kendall_w(rank_arr):
    """Kendall's coefficient of concordance W across per-patient rankings.

    Parameters
    ----------
    rank_arr : np.ndarray, shape (n_genes, n_patients)
        Column j is patient j's ranking of every gene (with ties averaged).

    Returns
    -------
    float in [0, 1] — 0 = no agreement, 1 = identical rankings; NaN if undefined.
    Tie-corrected per Siegel & Castellan.  We report the MAGNITUDE of W rather
    than its p-value because that p-value is meaningless at large n (any
    non-zero concordance becomes "significant").
    """
    n_items, m = rank_arr.shape
    if m < 2 or n_items < 2:
        return float('nan')
    Ri = rank_arr.sum(axis=1)
    S = float(((Ri - Ri.mean()) ** 2).sum())
    # tie correction term, summed over raters (patients)
    T = 0.0
    for j in range(m):
        _, counts = np.unique(rank_arr[:, j], return_counts=True)
        T += float(np.sum(counts ** 3 - counts))
    denom = m ** 2 * (n_items ** 3 - n_items) - m * T
    if denom <= 0:
        return float('nan')
    return float(12.0 * S / denom)


def _rra_scores(norm_rank_arr):
    """Robust Rank Aggregation scores (Kolde et al. 2012), pure Python.

    Parameters
    ----------
    norm_rank_arr : np.ndarray, shape (n_genes, n_patients)
        Normalized ranks in (0, 1]; smaller = better (ranked nearer the top).

    Returns
    -------
    (rho, score) : np.ndarray, np.ndarray
        For each gene, sort its normalized ranks ascending r_(1) <= ... <= r_(k).
        Under the null (uniform ranks) the i-th order statistic follows
        Beta(i, k-i+1); ``rho`` is the minimum of those Beta CDF probabilities
        over i.  ``score = min(rho * k, 1)`` applies a Bonferroni correction
        across the k order statistics.  Lower score = stronger, more consistent
        cross-cohort target.  Scores behave like p-values (used directly for BH).
    """
    n_genes, k = norm_rank_arr.shape
    if k == 0:
        return np.ones(n_genes), np.ones(n_genes)
    R = np.sort(norm_rank_arr, axis=1)                 # ascending order statistics
    i = np.arange(1, k + 1)
    P = beta.cdf(R, i, k - i + 1)                       # broadcast a=i, b=k-i+1
    rho = P.min(axis=1)
    score = np.minimum(rho * k, 1.0)
    return rho, score


def _robust_scale_per_patient(score_df):
    """Robustly-centered per-patient z-score of the gene scores (column-wise).

    Raw MOMA recovery scores are unbounded and dominated by a few large-flux
    reactions, so a patient's scores can reach ~1e11 while the bulk of genes sit
    in a near-tied baseline cluster. Testing such values against 0 is meaningless
    (0 is not the baseline) and saturates the rank-biserial/Wilcoxon statistics.

    We therefore center each patient on its MEDIAN knockout score — the
    per-patient control-centroid baseline, robust to the heavy upper tail — and
    standardize by the patient's standard deviation. Centering re-expresses each
    score as the SIGNED recovery relative to that baseline (positive = recovers
    more toward Control than the typical knockout), which de-saturates the
    sign/rank-based tests; the SD divisor (inflated by the few large outliers)
    keeps even the strongest targets on a modest scale rather than re-inflating
    them the way a tiny IQR would. Falls back to MAD, then 1.0, when the SD is 0;
    NaNs (genes absent from a patient) are preserved. Monotonic within each
    patient, so rank-based statistics (RRA, recurrence, Kendall's W) are
    unaffected — only median_score and the mixed-model estimates change scale.
    """
    scaled = score_df.copy()
    for col in scaled.columns:
        x = scaled[col]
        med = x.median()
        spread = x.std()                                        # SD (outlier-aware)
        if not np.isfinite(spread) or spread == 0:
            mad = (x - med).abs().median()
            spread = mad * 1.4826 if (np.isfinite(mad) and mad > 0) else 1.0
        scaled[col] = (x - med) / spread
    return scaled


def _normalized_rank_matrix(score_df):
    """Build per-patient rank and normalized-rank matrices from a score table.

    Higher score = better target (rank 1).  Genes absent from a patient's MOMA
    result are tied at the bottom (the average of the unused ranks), so every
    patient ranks the full gene union.  Returns (rank_df, norm_arr) where
    ``norm_arr`` holds normalized ranks in (0, 1].
    """
    G = score_df.shape[0]
    rank_df = score_df.rank(ascending=False, method='average')
    for col in rank_df.columns:
        m = int(rank_df[col].notna().sum())
        fill = (m + 1 + G) / 2.0 if m < G else np.nan
        rank_df[col] = rank_df[col].fillna(fill)
    norm_arr = (rank_df.values / G)
    return rank_df, norm_arr


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def rank_targets_across_cohort(scores_df, output_dir, effect_size_min=0.1,
                               fdr_cutoff=0.05, n_boot=1000, n_perm=1000,
                               seed=42, top_k=10, mixedlm_max_genes=300,
                               method='auto', switch_n=40):
    """Aggregate per-patient MOMA gene scores into a cohort-level ranking.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Full per-gene-per-patient table with at least columns
        ``Gene``, ``Sample``, ``score`` (and optionally ``Symbol``).  ``score``
        is the per-patient recovery metric (higher = better target).
    output_dir : str or None
        Directory to write ``cohort_target_ranking.csv``; if None, nothing is
        written (the DataFrame is still returned).
    effect_size_min : float
        Minimum effect size (recurrence proportion) required for significance.
    fdr_cutoff : float
        BH-FDR threshold for significance.
    n_boot, n_perm : int
        Bootstrap resamples (robustness) and label permutations (sensitivity).
    seed : int
        Master RNG seed — all randomness is derived from it (reproducible).
    top_k : int
        Top-k cutoff for recurrence and bootstrap selection frequency.
    mixedlm_max_genes : int
        Only fit the mixed model when the gene count does not exceed this — a
        per-gene fixed-effect model is degenerate and prohibitively slow for
        thousands of genes, so above the cap we fall back to RRA-only
        (mixedlm_estimate = NaN).
    method : {'auto', 'rra', 'mixedlm'}
        Inference strategy. RRA + permutation are the always-on backbone in every
        mode; this flag only controls whether the linear mixed-effects model is
        added.
          * 'auto'    — size-adaptive (default): add the mixed model only when it
            converges AND n_patients >= ``switch_n``; otherwise RRA only.
          * 'rra'     — force RRA only; never fit the mixed model.
          * 'mixedlm' — force the mixed model regardless of cohort size, still
            falling back to RRA gracefully if it fails to converge.
    switch_n : int
        Patient-count threshold at which 'auto' adds the mixed model (default 40).

    Returns
    -------
    pd.DataFrame, one row per gene, sorted ascending by ``rra_score`` (best
    first).  A gene is ``significant`` iff ``rra_fdr < fdr_cutoff`` AND
    ``recurrence_prop >= effect_size_min`` — FDR plus a mandatory effect size.
    """
    if method not in ('auto', 'rra', 'mixedlm'):
        raise ValueError("method must be one of 'auto', 'rra', 'mixedlm'")

    rng = np.random.default_rng(seed)

    required = {'Gene', 'Sample', 'score'}
    missing = required - set(scores_df.columns)
    if missing:
        raise ValueError('scores_df missing required columns: %s' % sorted(missing))

    df = scores_df.copy()
    df['Gene'] = df['Gene'].astype(str)
    df['Sample'] = df['Sample'].astype(str)
    df = df.dropna(subset=['score'])

    # genes x patients score matrix (mean if a gene appears twice for a patient)
    score_df = df.pivot_table(index='Gene', columns='Sample',
                              values='score', aggfunc='mean')

    # Per-patient robust scaling: tame unbounded fluxes and re-express each score
    # as the SIGNED recovery relative to the patient's control-centroid baseline,
    # before any aggregation or testing. Ranks (hence RRA/recurrence) are
    # preserved; the Wilcoxon/rank-biserial test and median_score are de-saturated
    # and put on a modest scale.
    score_df = _robust_scale_per_patient(score_df)
    # Scaled long form (used by the mixed model so its estimates are comparable).
    mixedlm_df = (score_df.reset_index()
                          .melt(id_vars='Gene', var_name='Sample', value_name='score')
                          .dropna(subset=['score']))

    genes = list(score_df.index)
    patients = list(score_df.columns)
    G = len(genes)
    n = len(patients)
    if G == 0 or n == 0:
        raise ValueError('scores_df contains no usable gene/patient scores.')
    k = min(int(top_k), G)

    # Symbol map (fall back to the gene id)
    if 'Symbol' in df.columns:
        sym = df.dropna(subset=['Symbol']).groupby('Gene')['Symbol'].first()
    else:
        sym = pd.Series(dtype=object)
    symbols = [str(sym.get(g, g)) for g in genes]

    # ---- (a) Agreement gate: Kendall's W ----------------------------------
    rank_df, norm_arr = _normalized_rank_matrix(score_df)
    W = _kendall_w(rank_df.values)
    low_agreement = bool(np.isnan(W) or W < 0.1)
    if low_agreement:
        logging.warning('Low cohort agreement (Kendall W = %.3f) — cross-patient '
                        'rankings are weakly concordant; interpret with care.', W)

    # ---- (b) Primary ranking: Robust Rank Aggregation ---------------------
    _, rra_score = _rra_scores(norm_arr)
    _, rra_fdr = fdrcorrection(rra_score)

    # ---- (c) Recurrence: top-k membership vs Binomial(n, k/G) -------------
    min_rank_df = score_df.rank(ascending=False, method='min')
    recurrence_count = min_rank_df.le(k).sum(axis=1).reindex(genes).fillna(0).astype(int).values
    recurrence_prop = recurrence_count / float(n)
    p0 = k / float(G)
    recurrence_p = binom.sf(recurrence_count - 1, n, p0)   # one-sided P(X >= count)
    _, recurrence_fdr = fdrcorrection(recurrence_p)

    # ---- (d) Effect consistency: Wilcoxon signed-rank of the per-patient
    #          robust-scaled signed recovery vs 0 (0 = control-centroid baseline).
    #          Scaling de-saturates rank-biserial and keeps median_score modest.
    wilcoxon_p = np.full(G, np.nan)
    rank_biserial = np.zeros(G)
    median_score = np.full(G, np.nan)
    for gi, g in enumerate(genes):
        vals = score_df.loc[g].values.astype(float)
        present = vals[~np.isnan(vals)]
        if present.size:
            median_score[gi] = float(np.median(present))
        d = np.where(np.isnan(vals), 0.0, vals)          # absent patient => 0 effect
        nz = d[d != 0]
        if nz.size < 1:
            continue
        ranks = rankdata(np.abs(nz))
        t_plus = ranks[nz > 0].sum()
        t_minus = ranks[nz < 0].sum()
        tot = t_plus + t_minus
        rank_biserial[gi] = (t_plus - t_minus) / tot if tot > 0 else 0.0
        try:
            wilcoxon_p[gi] = wilcoxon(nz)[1]
        except Exception:
            wilcoxon_p[gi] = np.nan
    wilcoxon_q = np.full(G, np.nan)
    qmask = ~np.isnan(wilcoxon_p)
    if qmask.any():
        _, q = fdrcorrection(wilcoxon_p[qmask])
        wilcoxon_q[qmask] = q

    # ---- (e) Robustness: bootstrap patients, top-k selection frequency ----
    boot_hits = np.zeros(G)
    if n_boot and n >= 2:
        for _ in range(int(n_boot)):
            idx = rng.integers(0, n, n)
            _, sc = _rra_scores(norm_arr[:, idx])
            boot_hits[np.argsort(sc)[:k]] += 1
        bootstrap_freq = boot_hits / float(n_boot)
    else:
        bootstrap_freq = np.full(G, np.nan)

    # ---- (f) Sensitivity: within-patient label-permutation null ----------
    if n_perm:
        null = np.empty(int(n_perm) * G)
        for p in range(int(n_perm)):
            # shuffle gene labels within each patient (independently per column),
            # vectorized: a random permutation of gene rows for every column at once
            perm_idx = np.argsort(rng.random(norm_arr.shape), axis=0)
            permuted = np.take_along_axis(norm_arr, perm_idx, axis=0)
            _, sc = _rra_scores(permuted)
            null[p * G:(p + 1) * G] = sc
        null_sorted = np.sort(null)
        # lower rra_score = stronger; empirical p = (1 + #{null <= obs}) / (1 + N)
        cnt = np.searchsorted(null_sorted, rra_score, side='right')
        perm_pvalue = (1.0 + cnt) / (1.0 + null.size)
    else:
        perm_pvalue = np.full(G, np.nan)

    # ---- (g) Size-adaptive switch: RRA/permutation are the always-on backbone;
    #          the mixed model is ADDED per `method` / `switch_n` (graceful fallback)
    mixedlm_estimate = np.full(G, np.nan)
    if method == 'rra':
        attempt_mixedlm = False
    elif method == 'mixedlm':
        attempt_mixedlm = True                      # forced regardless of cohort size
    else:  # 'auto'
        attempt_mixedlm = (n >= switch_n)           # add only for large enough cohorts
    if attempt_mixedlm:
        if G <= mixedlm_max_genes:
            mixedlm_estimate = _fit_mixedlm(mixedlm_df, genes)   # NaN on non-convergence
        else:
            logging.warning('Skipping mixed model: %d genes exceed mixedlm_max_genes=%d '
                            '(a per-gene fixed-effect model is degenerate at this size) '
                            '— using RRA only.', G, mixedlm_max_genes)
    elif method == 'auto':
        logging.info('Cohort target ranking: n=%d < switch_n=%d — RRA/permutation '
                     'only (mixed model not added).', n, switch_n)

    # ---- (i) GPR non-independence: identical-score gene groups -----------
    gpr_group = _gpr_groups(score_df, genes)

    # ---- (h) Significance = FDR AND effect size (recurrence proportion) ---
    significant = (rra_fdr < fdr_cutoff) & (recurrence_prop >= effect_size_min)

    out = pd.DataFrame({
        'Gene': genes,
        'Symbol': symbols,
        'gpr_group': gpr_group,
        'kendall_W': W,
        'low_cohort_agreement': low_agreement,
        'rra_score': rra_score,
        'rra_fdr': rra_fdr,
        'recurrence_count': recurrence_count,
        'recurrence_prop': recurrence_prop,
        'recurrence_fdr': recurrence_fdr,
        'wilcoxon_p': wilcoxon_p,
        'wilcoxon_q': wilcoxon_q,
        'rank_biserial': rank_biserial,
        'median_score': median_score,
        'bootstrap_freq': bootstrap_freq,
        'perm_pvalue': perm_pvalue,
        'mixedlm_estimate': mixedlm_estimate,
        'n_patients': n,
        'significant': significant,
    })
    out = out.sort_values('rra_score', kind='mergesort').reset_index(drop=True)

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        out.to_csv(os.path.join(output_dir, 'cohort_target_ranking.csv'), index=False)

    return out


def _fit_mixedlm(df, genes):
    """Per-gene fixed effects from a linear mixed model (random patient intercept).

    Fits ``score ~ gene`` with a per-patient random intercept via statsmodels
    MixedLM (built directly from design matrices to avoid the formula API).
    Returns an array of per-gene estimates relative to the reference gene
    (reference = 0.0), or all-NaN if the model is unavailable or fails to
    converge — the caller then relies on RRA/permutation alone.
    """
    G = len(genes)
    est = np.full(G, np.nan)
    try:
        from statsmodels.regression.mixed_linear_model import MixedLM
        long = df[['Gene', 'Sample', 'score']].dropna().copy()
        long['Gene'] = pd.Categorical(long['Gene'], categories=genes)
        dummies = pd.get_dummies(long['Gene'], drop_first=True).astype(float)
        exog = np.column_stack([np.ones(len(long)), dummies.values])
        model = MixedLM(long['score'].values, exog, groups=long['Sample'].values)
        res = model.fit()
        gene_idx = {g: i for i, g in enumerate(genes)}
        est[gene_idx[genes[0]]] = 0.0                     # reference level
        # fe_params: [Intercept, dummy_1, ...] aligned to dummies.columns
        fe = np.asarray(res.fe_params)
        for col, val in zip(dummies.columns, fe[1:]):
            if col in gene_idx:
                est[gene_idx[col]] = float(val)
        logging.info('Mixed model fitted for cohort target ranking (n>=30).')
    except Exception as e:
        logging.warning('Mixed model unavailable/failed (%s) — using RRA only.', e)
        est = np.full(G, np.nan)
    return est


def _gpr_groups(score_df, genes):
    """Tag genes with identical score vectors (GPR-collinear) into shared groups.

    Genes whose per-patient score vector is identical across the whole cohort are
    perturbed by the same reaction(s) and provide non-independent evidence.  Each
    multi-member identical-score group gets a label ``GPR_<k>``; singletons get an
    empty string.
    """
    G = len(genes)
    gpr_group = np.array([''] * G, dtype=object)
    keymap = {}
    for gi, g in enumerate(genes):
        row = score_df.loc[g].values
        key = tuple('NA' if (isinstance(v, float) and np.isnan(v)) else round(float(v), 9)
                    for v in row)
        keymap.setdefault(key, []).append(gi)
    gid = 0
    for members in keymap.values():
        if len(members) > 1:
            label = 'GPR_%d' % gid
            gid += 1
            for mi in members:
                gpr_group[mi] = label
    return gpr_group
