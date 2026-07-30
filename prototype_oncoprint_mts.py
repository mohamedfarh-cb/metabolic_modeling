#!/usr/bin/env python
"""PROTOTYPE — annotate the cohort target oncoprint with rMTA-style mTS excess.

STANDALONE. Not wired into run_pipeline; modifies nothing under tf_metabolism/.
All outputs are PROTO_-prefixed so they cannot overwrite pipeline figures.

WHY THIS EXISTS
---------------
Joining `cohort_target_ranking.csv` against `viz/prototype_rmta_mTS.csv` shows the
oncoprint's significant gene set and the rMTA transformation score disagree, and
disagree ASYMMETRICALLY:

    CD-ileum : 21 significant,  12 BELOW the cohort's own mTS background median
    UC-rectum: 40 significant,   2 below  (30 above the 98th percentile)

The worst case is CD glycolysis: ENO1/ENO2/ENO3 carry rra_fdr = 2.0e-89 and
recurrence_prop = 1.00 (top-10 lever in all 36 patients) yet sit at the 1.05th
percentile of their cohort's all-gene mTS background — i.e. knocking them out moves
CD patients toward control LESS than 99% of all genes in the model.

Mechanism: `recovery_proj_pca` is a signed projection but is NOT normalised by total
flux change, so perturbation magnitude still leaks into the score (the open [P0]
residual). mTS (rMTA Eq. 8) divides by total change, so magnitude cannot win it.

Nothing on the current oncoprint reveals any of this. This script makes it visible.

WHAT IT DOES / DOES NOT DO
--------------------------
DOES  : renders the EXISTING oncoprint (identical rows, identical cells) plus a
        far-right diverging bar of mTS excess over the cohort background, and a
        per-cohort all-gene concordance scatter of pipeline rank vs mTS excess.
DOES NOT: filter, re-rank, or re-select genes by mTS. Row selection is byte-identical
        to utils.visualize_target_oncoprint. An mTS-FILTERED oncoprint was considered
        and rejected (membership rule and cell values would encode contradicting
        metrics, and CD would collapse to 3 genes); an mTS-RANKED oncoprint (RRA over
        per-patient TS ranks) is the coherent alternative and remains deferred.
DOES NOT: constitute a significance test. mTS here is median-only, no bTS/wTS.
        [P0] stays open; nothing here validates a target.

Design: docs/superpowers/specs/2026-07-30-oncoprint-mts-annotation-design.md

Usage:
    python prototype_oncoprint_mts.py                       # both Human2 runs
    python prototype_oncoprint_mts.py --run CD=output_.../   # explicit
"""
from __future__ import print_function

import argparse
import copy as _copy
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib.lines import Line2D

# Style comes from the pipeline so the prototype is visually comparable to the real
# figure. Imported read-only — the same thing model_comparison.py does (CLAUDE.md 3e).
# Deliberately NOT wrapped in try/except: a silent fallback to default matplotlib style
# would make the prototype non-comparable to the figure it is meant to annotate.
from tf_metabolism.utils import set_publication_style


# Default targets: the two Human2 runs. The two Recon2M.2 runs have no
# viz/prototype_rmta_mTS.csv and are deliberately out of scope.
DEFAULT_RUNS = [
    ('CD_ileum', 'output_CD_ileum_Human2_tINIT_Human2_tINIT'),
    ('UC_rectum', 'output_UC_rectum_Human2_tINIT_Human2_tINIT'),
]

# Okabe-Ito. #0072B2 already means "KO" on the per-target cohort figure (3i).
POS_COLOR = '#0072B2'      # mTS excess above cohort background
NEG_COLOR = '#D55E00'      # below background — the finding this figure exists to show
MISSING_COLOR = '#BBBBBB'  # no mTS value (must not read as "at background")

# Fraction of the cohort's patients a gene must appear in to enter the background.
# REQUIRED: without it, genes present in only 1-2 patients' models top the ranking
# (CD: DCT/GNMT/SARDH; UC: ADSS2/TNRC6B/SLC17A7) — WORKLOG 2026-07-29c.
MIN_PATIENT_FRAC = 0.5


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def load_cohort(run_dir):
    """Join the cohort ranking with the mTS table; return (joined, filtered, background).

    `joined`   : every gene in cohort_target_ranking.csv, with mTS columns merged in
                 (NaN where the gene is absent from the mTS table).
    `filtered` : `joined` restricted to genes passing the min-patient filter — this is
                 the BACKGROUND population; every percentile is taken against it.
    `background`: median mTS_median over `filtered`.
    """
    ctr_csv = os.path.join(run_dir, 'cohort_target_ranking.csv')
    full_csv = os.path.join(run_dir, 'viz', 'moma_gene_scores_full.csv')
    mts_csv = os.path.join(run_dir, 'viz', 'prototype_rmta_mTS.csv')

    for f in (ctr_csv, full_csv):
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            raise SystemExit('FATAL missing/empty pipeline input: %s' % f)
    if not os.path.exists(mts_csv) or os.path.getsize(mts_csv) == 0:
        # Fail loud, with the fix. Rendering an unannotated figure here would silently
        # produce something indistinguishable from the real oncoprint.
        raise SystemExit(
            'FATAL no mTS table: %s\n'
            '      Generate it first:  python prototype_rmta_ts.py %s\n'
            '      (~10 min/cohort; reads the full stored MOMA matrices.)'
            % (mts_csv, run_dir))

    ctr = pd.read_csv(ctr_csv)
    full = pd.read_csv(full_csv)
    mts = pd.read_csv(mts_csv)
    if ctr.empty or full.empty or mts.empty:
        raise SystemExit('FATAL empty input table in %s' % run_dir)

    ctr['Gene'] = ctr['Gene'].astype(str)
    full['Gene'] = full['Gene'].astype(str)
    mts['Gene'] = mts['Gene'].astype(str)

    keep = mts[['Gene', 'mTS_median', 'n_patients']].rename(
        columns={'n_patients': 'n_patients_mTS'})
    joined = ctr.merge(keep, on='Gene', how='left')

    n_max = joined['n_patients_mTS'].max()
    passes = joined['n_patients_mTS'] >= MIN_PATIENT_FRAC * n_max
    filtered = joined[passes].dropna(subset=['mTS_median']).copy()
    background = float(filtered['mTS_median'].median())

    joined['background'] = background
    joined['excess'] = joined['mTS_median'] - background
    filtered['background'] = background
    filtered['excess'] = filtered['mTS_median'] - background

    # Percentile against the background population, not against the displayed genes.
    ref = np.sort(filtered['mTS_median'].values)

    def _pctile(v):
        if not np.isfinite(v):
            return np.nan
        return 100.0 * np.searchsorted(ref, v) / float(len(ref))

    joined['excess_pctile'] = joined['mTS_median'].apply(_pctile)
    filtered['excess_pctile'] = filtered['mTS_median'].apply(_pctile)

    # Genes the min-patient filter removed but that the ranking still calls significant:
    # they get a blank slot on the figure, so their excess must not be used.
    dropped = joined[(~passes) & joined['significant'].astype(str)
                     .str.lower().isin(['true', '1', '1.0'])]
    if len(dropped):
        joined.loc[~passes, ['excess', 'excess_pctile']] = np.nan

    return joined, filtered, background, full


def select_rows(ctr, max_genes=40, fallback_n=20):
    """Row selection + ordering, byte-identical to utils.visualize_target_oncoprint.

    Kept verbatim on purpose: the prototype must annotate exactly the rows production
    shows. Any divergence here would make the annotation answer a different question.
    """
    sig_mask = ctr['significant'].astype(str).str.lower().isin(['true', '1', '1.0'])
    fallback = False
    if sig_mask.any():
        sel = ctr[sig_mask].copy()
        if len(sel) > max_genes:
            sel = sel.sort_values('rra_score', kind='mergesort').head(max_genes)
    else:
        fallback = True
        sel = ctr.sort_values('rra_score', kind='mergesort').head(fallback_n).copy()

    sel['gpr_group'] = sel['gpr_group'].fillna('').astype(str)
    sel['_grp'] = [g.strip() if g.strip() else '__solo_%s' % gene
                   for g, gene in zip(sel['gpr_group'], sel['Gene'])]
    sel = sel.sort_values(['recurrence_prop', 'rra_score'],
                          ascending=[False, True], kind='mergesort').reset_index(drop=True)
    sel['_pos'] = np.arange(len(sel))
    sel['_grp_first'] = sel.groupby('_grp')['_pos'].transform('min')
    sel = sel.sort_values(['_grp_first', '_pos'], kind='mergesort').reset_index(drop=True)
    return sel, fallback


# ---------------------------------------------------------------------------
# Figure 1 — oncoprint + mTS excess bar
# ---------------------------------------------------------------------------
def render_oncoprint(cohort, run_dir, sel, full, fallback, background,
                     top_k=10, fallback_n=20):
    """Copy of utils.visualize_target_oncoprint with a far-right mTS-excess panel.

    COPY, NOT IMPORT: the production function saves its figure and returns nothing, and
    the new panel requires editing the fixed add_axes fractions from inside the layout
    block. This copy will drift from utils.py — acceptable for a go/no-go prototype, and
    it guarantees zero risk to the pipeline (same trade-off prototype_ko_cohort_calibrated.R
    made against moma_cohort_viz.R).
    """
    out_viz = os.path.join(run_dir, 'viz')
    k = max(int(top_k), 1)

    row_genes = list(sel['Gene'])

    def _label(row):
        s = row.get('Symbol')
        return str(s) if isinstance(s, str) and s and s.lower() != 'nan' else str(row['Gene'])

    row_labels = [_label(r) for _, r in sel.iterrows()]
    row_recur = sel['recurrence_prop'].astype(float).values
    row_deflux = (sel['de_flux'].astype(str).str.lower().isin(['true', '1', '1.0']).values
                  if 'de_flux' in sel.columns else np.zeros(len(sel), dtype=bool))
    row_organelle = (sel['organelle'].astype(str).values if 'organelle' in sel.columns
                     else np.array([''] * len(sel), dtype=object))
    row_excess = sel['excess'].astype(float).values          # NEW

    ORG_PAL = {'Cytosol': '#555555', 'Mitochondria': '#D55E00', 'Peroxisome': '#009E73',
               'ER': '#0072B2', 'Golgi': '#CC79A7', 'Nucleus': '#E69F00',
               'Lysosome': '#56B4E9', 'Inner mito': '#D55E00', 'Extracellular': '#999999'}

    def _org_color(o):
        return ORG_PAL.get(str(o), '#555555')

    # ---- Rank matrix (genes x patients) ---------------------------------------
    fsel = full[full['Gene'].isin(set(row_genes))]
    rankmat = fsel.pivot_table(index='Gene', columns='Sample', values='rank', aggfunc='min')
    patients_all = list(full['Sample'].astype(str).unique())
    rankmat = rankmat.reindex(index=row_genes, columns=patients_all)

    M = rankmat.values.astype(float)
    M[~np.isfinite(M)] = np.nan
    M[M > k] = np.nan

    col_hits_all = np.isfinite(M).sum(axis=0)
    order = np.argsort(-col_hits_all, kind='stable')
    patients = [patients_all[i] for i in order]
    M = M[:, order]
    col_hits = np.isfinite(M).sum(axis=0)
    n, m = M.shape

    masked = np.ma.masked_invalid(M)
    cmap = _copy.copy(cm.get_cmap('viridis'))
    cmap.set_bad('#f0f0f0')
    norm = Normalize(vmin=1, vmax=max(k, 2))

    # ---- Layout: ax_main narrowed, ax_right narrowed, ax_mts inserted ---------
    set_publication_style()
    fig_w = max(9.0, 0.24 * m + 5.8)             # +1.0 floor vs production for the new column
    fig_h = max(4.5, 0.32 * n + 2.6)
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')
    fig.set_constrained_layout(False)
    ax_left = fig.add_axes([0.015, 0.085, 0.150, 0.70])
    ax_main = fig.add_axes([0.175, 0.085, 0.545, 0.70])
    ax_top = fig.add_axes([0.175, 0.800, 0.545, 0.130], sharex=ax_main)
    ax_right = fig.add_axes([0.725, 0.085, 0.105, 0.70], sharey=ax_main)
    ax_mts = fig.add_axes([0.845, 0.085, 0.135, 0.70], sharey=ax_main)   # NEW
    ax_cbar = fig.add_axes([0.725, 0.840, 0.105, 0.028])

    im = ax_main.imshow(masked, aspect='auto', cmap=cmap, norm=norm,
                        origin='upper', interpolation='nearest')
    ax_main.set_xticks(np.arange(m))
    ax_main.set_xticklabels(patients, rotation=90, fontsize=max(6, 9 - m // 12))
    ax_main.set_yticks([])
    ax_main.set_xlim(-0.5, m - 0.5)
    ax_main.set_ylim(n - 0.5, -0.5)
    ax_main.set_xticks(np.arange(-0.5, m, 1), minor=True)
    ax_main.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax_main.grid(which='minor', color='white', linewidth=0.6)
    ax_main.tick_params(which='minor', length=0)

    # Left: gene labels + GPR brackets
    ax_left.set_xlim(0, 1)
    ax_left.set_ylim(n - 0.5, -0.5)
    ax_left.axis('off')
    for i, (lab, de, org) in enumerate(zip(row_labels, row_deflux, row_organelle)):
        txt = lab + (' *' if de else '')
        ax_left.text(0.96, i, txt, ha='right', va='center',
                     fontsize=max(6, 9 - n // 18),
                     fontweight='bold' if de else 'normal',
                     color=_org_color(org), fontstyle='italic')

    def _group_name(symbols):
        pref = os.path.commonprefix(symbols)
        pref = ''.join(ch for ch in pref if ch.isalpha())
        return pref if len(pref) >= 3 else '/'.join(symbols)

    def _short_pathway(s):
        s = str(s).split(' / ')[0].strip()
        return '' if s.lower() in ('', 'nan') else (s if len(s) <= 22 else s[:21] + u'…')

    grp_sub_enr = (sel['subsystem'].values if 'subsystem' in sel.columns
                   else np.array([''] * n, dtype=object))
    grp_sub_frq = (sel['subsystem_frequent'].values if 'subsystem_frequent' in sel.columns
                   else grp_sub_enr)
    grp_arr = sel['_grp'].values
    i = 0
    while i < n:
        j = i
        while j + 1 < n and grp_arr[j + 1] == grp_arr[i]:
            j += 1
        if j > i and not str(grp_arr[i]).startswith('__solo_'):
            xb = 0.17
            ax_left.plot([xb, xb], [i - 0.42, j + 0.42], color='#444444', lw=1.1, clip_on=False)
            ax_left.plot([xb, xb + 0.045], [i - 0.42, i - 0.42], color='#444444',
                         lw=1.1, clip_on=False)
            ax_left.plot([xb, xb + 0.045], [j + 0.42, j + 0.42], color='#444444',
                         lw=1.1, clip_on=False)
            enr = _short_pathway(grp_sub_enr[i])
            frq = _short_pathway(grp_sub_frq[i])
            lines = []
            if enr:
                lines.append(u'– %s –' % enr)
            if frq:
                lines.append('(%s)' % frq)
            glabel = '\n'.join(lines) if lines else _group_name(row_labels[i:j + 1])
            ax_left.text(0.0, (i + j) / 2.0, glabel, rotation=90, ha='left', va='center',
                         fontsize=6, color='#333333', linespacing=1.15)
        i = j + 1

    # Top: per-patient hit count
    ax_top.bar(np.arange(m), col_hits, color='#4C72B0', width=0.9, edgecolor='none')
    ax_top.set_ylabel('targets\n(top-%d)' % k, fontsize=7)
    ax_top.tick_params(axis='x', labelbottom=False, length=0)
    ax_top.tick_params(axis='y', labelsize=6)
    ax_top.set_ylim(0, max(col_hits.max(), 1) * 1.15)
    for s in ('top', 'right'):
        ax_top.spines[s].set_visible(False)

    # Right: per-gene recurrence proportion
    ax_right.barh(np.arange(n), row_recur, color='#DD8452', height=0.8, edgecolor='none')
    ax_right.set_ylim(n - 0.5, -0.5)
    ax_right.set_xlim(0, 1)
    ax_right.set_xlabel('recurrence', fontsize=7)
    ax_right.set_yticks([])
    ax_right.tick_params(axis='x', labelsize=6)
    ax_right.set_xticks([0, 0.5, 1.0])
    for s in ('top', 'right'):
        ax_right.spines[s].set_visible(False)

    # ---- NEW: mTS excess over cohort background -------------------------------
    # Adjacent to recurrence on purpose: the two competing effect sizes sit side by
    # side, so a disagreement is one horizontal scan per row.
    finite = np.isfinite(row_excess)
    lim = float(np.nanmax(np.abs(row_excess))) * 1.18 if finite.any() else 1.0
    idx = np.arange(n)
    bar_colors = [POS_COLOR if v >= 0 else NEG_COLOR for v in row_excess[finite]]
    ax_mts.barh(idx[finite], row_excess[finite], color=bar_colors,
                height=0.8, edgecolor='none')
    for i in idx[~finite]:
        # Blank slot, not a zero bar: no measurement must never read as "at background".
        ax_mts.plot([0], [i], marker='|', markersize=7, color=MISSING_COLOR,
                    markeredgewidth=1.2)
    ax_mts.axvline(0, color='#444444', lw=0.8)
    ax_mts.set_ylim(n - 0.5, -0.5)
    ax_mts.set_xlim(-lim, lim)
    ax_mts.set_xlabel('mTS excess\n(cohort bg %+.3f)' % background, fontsize=7)
    ax_mts.set_yticks([])
    ax_mts.set_xticks([-lim, 0, lim])
    ax_mts.set_xticklabels(['%+.2f' % -lim, '0', '%+.2f' % lim])
    ax_mts.tick_params(axis='x', labelsize=6)
    # The outer ticks sit ON the axis edges, so centre-aligned labels overhang and the
    # left one collides with the recurrence axis's "1.0" (seen at UC's lim=0.11). Pull
    # the two outer labels inward; the middle "0" stays centred on the zero line.
    _xt = ax_mts.get_xticklabels()
    _xt[0].set_ha('left')
    _xt[-1].set_ha('right')
    for s in ('top', 'right'):
        ax_mts.spines[s].set_visible(False)

    # Rank colorbar
    cb = fig.colorbar(im, cax=ax_cbar, orientation='horizontal')
    cb.set_label('within-patient rank (1 = top)', fontsize=7)
    cb.ax.tick_params(labelsize=6)
    cb.set_ticks([1, max(k, 2)])

    # Titles / legends
    fig.suptitle('PROTOTYPE — cohort target oncoprint with mTS calibration (%s)'
                 % cohort.replace('_', '-'), x=0.5, y=0.985, fontsize=13,
                 fontweight='medium')
    n_de = int(np.count_nonzero(row_deflux))
    n_below = int(np.count_nonzero(row_excess[finite] < 0))
    sub = ('no significant targets — showing top %d by RRA' % fallback_n if fallback
           else ('%d significant target genes (all FDR + effect-size significant)  '
                 u'•  * = gene is DE-flux included (%d/%d)  '
                 u'•  %d/%d BELOW cohort mTS background'
                 % (n, n_de, n, n_below, int(finite.sum()))))
    fig.text(0.175, 0.945, sub, fontsize=8.5, color='#333333')

    conv_handles = [Line2D([0], [0], marker='|', color='#444444', lw=0,
                           markersize=11, markeredgewidth=1.4, label='GPR-coupled group'),
                    Line2D([0], [0], marker='$*$', color='black', lw=0,
                           markersize=10, label='*  DE-flux gene'),
                    Line2D([0], [0], color='none', label=u'– x –   most-enriched pathway'),
                    Line2D([0], [0], color='none', label='( x )    most-frequent subsystem'),
                    Line2D([0], [0], marker='s', color='none', markerfacecolor=POS_COLOR,
                           markersize=11, label='mTS above cohort background'),
                    Line2D([0], [0], marker='s', color='none', markerfacecolor=NEG_COLOR,
                           markersize=11, label='mTS BELOW cohort background')]
    leg1 = fig.legend(handles=conv_handles, loc='upper left', bbox_to_anchor=(0.015, -0.015),
                      frameon=False, fontsize=9, handletextpad=0.6, labelspacing=0.6,
                      title='plot key', title_fontsize=9)
    leg1.get_title().set_fontweight('bold')
    fig.add_artist(leg1)

    seen_org = [o for o in dict.fromkeys(row_organelle) if str(o) and str(o) != 'nan']
    org_handles = [Line2D([0], [0], marker='s', color='none', markerfacecolor=_org_color(o),
                          markersize=11, label=o)
                   for o in sorted(seen_org, key=lambda x: list(ORG_PAL).index(x)
                                   if x in ORG_PAL else 99)]
    if org_handles:
        leg2 = fig.legend(handles=org_handles, loc='upper left', bbox_to_anchor=(0.42, -0.015),
                          frameon=False, fontsize=9, handletextpad=0.5, labelspacing=0.6,
                          ncol=2, columnspacing=1.4,
                          title='gene-label colour = organelle', title_fontsize=9)
        leg2.get_title().set_fontweight('bold')

    n_missing = int((~finite).sum())
    caveat = ('mTS = rMTA (Valcarcel 2019) Eq. 8, median over patients; no bTS/wTS,\n'
              'not a significance test. [P0] competitive null still open —\n'
              'hypothesis-only, never a validated target.')
    if n_missing:
        caveat += '\n%d shown gene(s) lack an mTS value (grey tick, not zero).' % n_missing
    fig.text(0.72, -0.02, caveat, fontsize=7, color='#666666', va='top', ha='left')

    base = os.path.join(out_viz, 'PROTO_target_oncoprint_mts_%s' % cohort)
    fig.savefig(base + '.png', format='png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(base + '.svg', format='svg', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return base + '.png', n_below, int(finite.sum()), n_missing


# ---------------------------------------------------------------------------
# Figure 2 — all-gene concordance scatter
# ---------------------------------------------------------------------------
def render_scatter(cohort, run_dir, filtered, background, rho_rank):
    """Pipeline rank (x, log) vs mTS excess (y) over ALL background genes.

    Answers what the annotation bar structurally cannot: whether the two metrics
    disagree only among the significant genes or across the whole model.
    """
    out_viz = os.path.join(run_dir, 'viz')
    d = filtered.copy()
    d['rra_rank'] = d['rra_score'].rank(method='min', ascending=True)
    sig = d[d['significant'].astype(str).str.lower().isin(['true', '1', '1.0'])] \
        .sort_values('excess', ascending=False)

    set_publication_style()
    fig, ax = plt.subplots(figsize=(7.6, 6.0), facecolor='white')
    ax.axhline(0, color='#444444', lw=0.9, zorder=1)
    ax.scatter(d['rra_rank'], d['excess'], s=6, color='#CCCCCC', linewidths=0,
               zorder=2, label='all genes (n=%d)' % len(d))
    ax.scatter(sig['rra_rank'], sig['excess'], s=34,
               color=[POS_COLOR if v >= 0 else NEG_COLOR for v in sig['excess']],
               edgecolor='white', linewidths=0.6, zorder=4,
               label='oncoprint significant (n=%d)' % len(sig))

    # ---- Robust y-window ------------------------------------------------------
    # A handful of extreme background genes (CD reaches -0.4) otherwise stretch the
    # axis and squash every significant gene into the top 15% of the panel. Clip to a
    # window covering all significant genes plus the 1st-99th percentile of the
    # background, and SAY how many points fell outside — a silently cropped axis would
    # misrepresent the spread of the background this figure exists to characterise.
    span = np.concatenate([sig['excess'].values,
                           np.nanpercentile(d['excess'].values, [1, 99])])
    half = float(np.nanmax(np.abs(span))) * 1.30
    ax.set_ylim(-half, half)
    n_clipped = int((np.abs(d['excess']) > half).sum())

    # ---- Deterministic label placement ---------------------------------------
    # adjustText is NOT installed here, and a seeded jitter would break the repo's
    # determinism convention.
    #
    # KEY POINT: GPR-collinear genes share an IDENTICAL knockout effect, so they have an
    # IDENTICAL excess by construction (PKM/PKLR, ENO1/2/3, NME1/2, GAPDH/GAPDHS...).
    # Labelling them separately guarantees overlapping text no matter how the labels are
    # nudged. Merge co-located genes into ONE label, which also makes the collinearity
    # visible rather than hiding it behind stacked names.
    #
    # Labels then sit directly beside their point with a small offset in POINTS, so they
    # cannot escape the data limits (an earlier data-coordinate cascade pushed them past
    # the axes and over the title). Only the 5 highest and 5 lowest groups are labelled —
    # per-gene detail lives on the oncoprint bar and in the concordance CSV.
    def _sym(r):
        return (r['Symbol'] if isinstance(r['Symbol'], str) and r['Symbol'].lower() != 'nan'
                else r['Gene'])

    sig = sig.copy()
    sig['_key'] = sig['excess'].round(6)
    groups = []
    for key, g in sig.groupby('_key', sort=False):
        syms = [_sym(r) for _, r in g.iterrows()]
        # Collapse a shared-prefix family to "ENO1/2/3" rather than "ENO1/ENO2/ENO3".
        pref = os.path.commonprefix(syms)
        sufs = [s[len(pref):] for s in syms]
        # all(sufs): when one symbol IS the prefix its suffix is empty, which would emit a
        # trailing slash ("GAPDHS/" for GAPDHS+GAPDH). Fall back to the plain join there.
        if len(syms) > 1 and len(pref) >= 3 and all(sufs):
            label = pref + '/'.join(sufs)
        else:
            label = '/'.join(syms)
        if len(label) > 22:
            label = '%s +%d more' % (syms[0], len(syms) - 1)
        groups.append((float(key), float(np.median(g['rra_rank'].values)), label))
    groups.sort(key=lambda t: -t[0])

    n_lab = min(5, len(groups))
    labelled = groups[:n_lab] + groups[-n_lab:] if len(groups) > n_lab else groups
    seen_keys = set()
    for i, (excess_v, rank_v, label) in enumerate(labelled):
        if excess_v in seen_keys:
            continue
        seen_keys.add(excess_v)
        # Alternate the offset magnitude: the top groups sit at nearly the same y (UC's
        # five oxidase genes are all ~+0.085), so a constant offset stacks their text on
        # one line (ACOX3 over ACOX1, XDH over DAO).
        mag = 9 if i % 2 == 0 else 20
        ax.annotate(label, xy=(rank_v, excess_v),
                    xytext=(7, mag if excess_v >= 0 else -mag),
                    textcoords='offset points', fontsize=7, va='center', ha='left',
                    color='#333333', zorder=5)

    ax.set_xscale('log')
    ax.set_xlabel('pipeline rank by rra_score  (1 = strongest target)  [log]')
    ax.set_ylabel('mTS excess over cohort background')
    ax.set_title('PROTOTYPE — pipeline rank vs rMTA mTS excess (%s)'
                 % cohort.replace('_', '-'))
    note = (u'Spearman ρ(rra rank, −excess) = %+.3f\ncohort mTS background = %+.4f\n'
            u'labelled: %d highest / %d lowest groups (GPR-identical genes merged)'
            % (rho_rank, background, n_lab, n_lab))
    if n_clipped:
        note += '\n%d background gene(s) outside the y-window' % n_clipped
    # Below the axes as a caption, not inside them: whichever corner it occupies it
    # collides with labelled points in one cohort or the other (top-left buried UC's
    # NOS2, bottom-left buried CD's ENO group).
    fig.text(0.02, -0.005, note, fontsize=8, va='top', ha='left', color='#333333')
    ax.legend(loc='upper right', frameon=False, fontsize=8, markerscale=1.2)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)

    base = os.path.join(out_viz, 'PROTO_mts_concordance_scatter_%s' % cohort)
    fig.savefig(base + '.png', format='png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(base + '.svg', format='svg', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return base + '.png'


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', action='append', metavar='LABEL=DIR',
                    help='cohort to render (repeatable); default = the two Human2 runs')
    ap.add_argument('--out-csv', default='PROTO_oncoprint_mts_concordance.csv')
    ap.add_argument('--top_k', type=int, default=10, help='oncoprint cell cutoff (default 10)')
    ap.add_argument('--max_genes', type=int, default=40)
    ap.add_argument('--fallback_n', type=int, default=20)
    args = ap.parse_args()

    runs = DEFAULT_RUNS
    if args.run:
        runs = []
        for spec in args.run:
            if '=' not in spec:
                raise SystemExit('FATAL --run needs LABEL=DIR, got: %s' % spec)
            label, _, path = spec.partition('=')
            runs.append((label.strip(), path.strip().rstrip('/')))

    rows = []
    for cohort, run_dir in runs:
        if not os.path.isdir(run_dir):
            raise SystemExit('FATAL not a directory: %s' % run_dir)
        print('=== %s  (%s)' % (cohort, run_dir))

        joined, filtered, background, full = load_cohort(run_dir)
        sel, fallback = select_rows(joined, args.max_genes, args.fallback_n)

        png1, n_below, n_scored, n_missing = render_oncoprint(
            cohort, run_dir, sel, full, fallback, background,
            top_k=args.top_k, fallback_n=args.fallback_n)

        # Two DISTINCT Spearman statistics — they order the cohorts oppositely, so
        # quoting the wrong one flips which cohort looks more concordant.
        rho_rank, p_rank = spearmanr(filtered['rra_score'].rank(),
                                     (-filtered['excess']).rank())
        rho_score, p_score = spearmanr(filtered['median_score'], filtered['excess'])

        png2 = render_scatter(cohort, run_dir, filtered, background, rho_rank)

        pct = filtered['excess_pctile']
        sigf = filtered[filtered['significant'].astype(str).str.lower()
                        .isin(['true', '1', '1.0'])]
        spct = sigf['excess_pctile']
        print('    background median mTS   : %+.4f  (over %d genes, >=%d%% patients)'
              % (background, len(filtered), int(MIN_PATIENT_FRAC * 100)))
        print('    significant genes       : %d' % len(sigf))
        print('    below background        : %d' % int((sigf['excess'] < 0).sum()))
        print('    below 10th pctile       : %d' % int((spct < 10).sum()))
        print('    above 98th pctile       : %d' % int((spct > 98).sum()))
        print('    rho_rank  (rra rank vs -excess) : %+.3f  (p=%.2g)' % (rho_rank, p_rank))
        print('    rho_score (median_score vs excess): %+.3f  (p=%.2g)' % (rho_score, p_score))
        if n_missing:
            print('    shown genes lacking mTS : %d' % n_missing)
        print('    wrote %s' % png1)
        print('    wrote %s' % png2)

        keep = ['Gene', 'Symbol', 'gpr_group', 'de_flux', 'rra_score', 'rra_fdr',
                'recurrence_prop', 'median_score', 'mTS_median', 'background',
                'excess', 'excess_pctile', 'n_patients_mTS']
        t = sel[[c for c in keep if c in sel.columns]].copy()
        t.insert(0, 'cohort', cohort)
        rows.append(t)

    out = pd.concat(rows, ignore_index=True)
    out = out.sort_values(['cohort', 'rra_score'], kind='mergesort')
    out.to_csv(args.out_csv, index=False)
    print('\nwrote %s  (%d rows)' % (args.out_csv, len(out)))


if __name__ == '__main__':
    sys.exit(main())
