#!/usr/bin/env python
"""PROTOTYPE — patient-to-control distance in flux space, before vs after a knockout.

STANDALONE. Not wired into run_pipeline; modifies nothing under tf_metabolism/.
Outputs are PROTO_-prefixed so they cannot overwrite pipeline figures.

WHY THIS EXISTS
---------------
The UMAP cannot answer "did this knockout help" -- it compresses the 244 (CD) / 771 (UC)
DE-flux reactions to 2 numbers, discarding the per-reaction directions the question is
about (measured: departure-from-drift correlates with mTS at +0.006 CD / +0.159 UC).

rMTA's mTS fixes direction but answers a different question. mTS is defined on the CHANGE
(KO - WT) and is NORMALISED by total change, so it reports what FRACTION of the movement
helped -- never HOW MUCH distance was closed. A knockout can score mTS +0.9 while moving
almost nothing. This restores the missing magnitude by measuring STATES:

    distance(v) = sum_i  |v_i - c_i| / gap_i

  i      DE-flux reactions with Status in {UP, DOWN}; Reversed is DROPPED (no defined
         disease direction) and the count is reported on the figure.
  c_i    control mean          (Differential_fluxes.csv `Condition1 mean`)
  gap_i  |disease mean - control mean|, i.e. that reaction's own disease-vs-control
         separation. Distance is therefore read in DISEASE-GAPS: 0 = at the control
         mean, 1 = a typical diseased distance away.

distance(control mean) == 0 by construction and distance >= 0 always -- both asserted.

TWO DESIGN CHOICES THAT WERE WRONG IN THE FIRST VERSION (kept as options, not defaults)
--------------------------------------------------------------------------------------
1. `--deviation disease` rectifies to max(0, disease-direction deviation), counting only
   deviation on the disease side. It was the first default, on the reasoning that
   overshooting past control should not count as extra disease. That is backwards: the
   rectifier credits overshoot as a PERFECT correction, so knockouts that blow past
   control score best. Measured on UC patient C3003: a typical knockout overshoots 59
   reactions vs ACOX1's 49, and the rectified metric consequently ranked ACOX1 BELOW
   typical while the true distance ranked it far above. Default is now `abs`, which
   penalises overshoot like any distance.

2. `--scale sd` divides by the control SD. 52 UC reactions have control SD < 1e-6, and
   after scaling a SINGLE reaction accounted for 99.5% of the total distance -- the metric
   measured one reaction. Default is now `gap`, where the top reaction contributes
   6.9% (UC) / 9.9% (CD). `--scale raw` uses flux units and is dominated by high-flux
   reactions instead. Concentration is reported at run time either way.

CALIBRATION IS BUILT IN, NOT OPTIONAL
-------------------------------------
An uncalibrated count misleads (CD: 79% of ALL genes move toward control). So every panel
also plots the TYPICAL knockout: the median distance over ALL ~2,200 genes in that same
patient. A gene must beat its own grey point, not merely fall below the diagonal. This is
nearly free -- each patient's MOMA file is a reactions x genes matrix, so all genes are
scored in one vectorised pass.

Nothing is re-simulated; this reads stored MOMA output.

Usage:
    python prototype_control_distance.py
    python prototype_control_distance.py --genes ACOX1,CPT2 --scale raw
"""
from __future__ import print_function

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from tf_metabolism.utils import set_publication_style


DEFAULT_RUNS = [
    ('CD_ileum', 'output_CD_ileum_Human2_tINIT_Human2_tINIT'),
    ('UC_rectum', 'output_UC_rectum_Human2_tINIT_Human2_tINIT'),
]
DEFAULT_GENES = ['SLC7A2', 'ENO1', 'ACOX1', 'CPT2']

GENE_COLOR = '#0072B2'
TYPICAL_COLOR = '#BBBBBB'
WORSE_COLOR = '#D55E00'


def disease_axis(run_dir, scale_mode):
    """Per DE reaction: control mean, disease direction, per-reaction scale, n dropped."""
    df = pd.read_csv(os.path.join(run_dir, 'Differential_fluxes.csv'), index_col=0)
    status = df['Status'].astype(str).str.upper()
    keep = status.isin(['UP', 'DOWN'])
    n_dropped = int((~keep).sum())
    df, status = df[keep], status[keep]

    c = df['Condition1 mean'].astype(float)
    # +1 where disease sits ABOVE control. Negation of prototype_rmta_ts.desired_direction(),
    # which returns the direction of TRAVEL toward control. Used only by --deviation disease.
    s = pd.Series(np.where(status == 'UP', 1.0, -1.0), index=df.index)

    if scale_mode == 'gap':
        scale = (df['Condition2 mean'].astype(float) - c).abs()
    elif scale_mode == 'sd':
        sd = df['Condition1 std'].astype(float).abs()
        scale = sd.where(sd > 0, c.abs().where(c.abs() > 0, 1.0))
    else:
        scale = pd.Series(1.0, index=df.index)
    # A zero scale would divide by zero and swamp everything. For DE reactions the gap is
    # non-zero by construction (verified: 0 zero-gaps in both cohorts); this is a guard.
    scale = scale.where(scale > 0, 1.0)
    return c, s, scale, n_dropped


def distance_matrix(V, c, s, scale, deviation):
    """Distance to control for every column of V (reactions x items)."""
    dev = V.values - c.values[:, None]
    if deviation == 'disease':
        dev = np.maximum(dev * s.values[:, None], 0.0)   # see docstring: biased, non-default
    else:
        dev = np.abs(dev)
    return (dev / scale.values[:, None]).sum(axis=0)


def score_cohort(cohort, run_dir, gene_ids, scale_mode, deviation):
    c, s, scale, n_dropped = disease_axis(run_dir, scale_mode)
    wt_all = pd.read_csv(os.path.join(run_dir, 'flux2.csv'), index_col=0).fillna(0.0)
    files = sorted(glob.glob(os.path.join(run_dir, 'targeting_results',
                                          'MOMA_target_results_*.csv')))
    if not files:
        raise SystemExit('FATAL no MOMA output in %s/targeting_results' % run_dir)

    rows, n_rxn, top_share = [], None, None
    for n, f in enumerate(files, 1):
        sample = os.path.basename(f).replace('MOMA_target_results_', '').replace('.csv', '')
        if sample not in wt_all.columns:
            print('    skip %s (no WT flux column)' % sample)
            continue
        ko = pd.read_csv(f, index_col=0)
        rxns = [r for r in c.index if r in ko.index and r in wt_all.index]
        if not rxns:
            print('    skip %s (no DE reactions in its MOMA output)' % sample)
            continue
        n_rxn = len(rxns)
        cc, ss, sk = c.loc[rxns], s.loc[rxns], scale.loc[rxns]

        wt_vec = wt_all.loc[rxns, sample].fillna(0.0).to_frame('WT')
        d_wt = float(distance_matrix(wt_vec, cc, ss, sk, deviation)[0])
        d_all = pd.Series(distance_matrix(ko.loc[rxns].astype(np.float64), cc, ss, sk,
                                          deviation), index=ko.columns)
        d_typ = float(np.nanmedian(d_all.values))

        if top_share is None:      # concentration diagnostic, first usable patient
            per = np.abs(wt_vec.values[:, 0] - cc.values) / sk.values
            top_share = 100.0 * per.max() / per.sum() if per.sum() > 0 else np.nan

        for gid, label in gene_ids:
            if gid not in d_all.index:
                continue
            d_ko = float(d_all[gid])
            # Signed so POSITIVE = moved AWAY from control. Every knockout measured so far
            # is positive (see the docstring finding), which is why this is not phrased as
            # "% closed" -- that label would print negative numbers for every gene.
            rows.append(dict(cohort=cohort, sample=sample, gene=label,
                             dist_wt=d_wt, dist_ko=d_ko, dist_typical_ko=d_typ,
                             pct_change=100.0 * (d_ko - d_wt) / d_wt if d_wt > 0 else np.nan,
                             pct_change_typical=100.0 * (d_typ - d_wt) / d_wt
                             if d_wt > 0 else np.nan))
        print('    [%d/%d] %s: %d reactions, WT %.1f, typical KO %.1f'
              % (n, len(files), sample, len(rxns), d_wt, d_typ))
        sys.stdout.flush()

    return pd.DataFrame(rows), n_rxn, n_dropped, top_share


def render(cohort, run_dir, res, scale_mode, deviation, n_rxn, n_dropped, labels):
    labels = [l for l in labels if l in set(res['gene'])]
    set_publication_style()
    fig, axes = plt.subplots(1, len(labels), figsize=(3.5 * len(labels) + 1.0, 4.6),
                             facecolor='white', squeeze=False)
    axes = axes[0]
    hi = max(res['dist_wt'].max(), res['dist_ko'].max(),
             res['dist_typical_ko'].max()) * 1.08
    lo = -0.02 * hi

    for ax, lab in zip(axes, labels):
        g = res[res['gene'] == lab]
        ax.plot([lo, hi], [lo, hi], color='#888888', lw=1.0, ls='--', zorder=1)
        ax.scatter(g['dist_wt'], g['dist_typical_ko'], s=22, color=TYPICAL_COLOR,
                   linewidths=0, zorder=2, label='typical KO (all-gene median)')
        worse = g['dist_ko'] > g['dist_wt']
        # Only label a category that HAS points. Every knockout measured so far moves the
        # patient away, so an unconditional "<gene> KO" swatch advertises a blue series
        # the panel does not contain.
        if (~worse).any():
            ax.scatter(g.loc[~worse, 'dist_wt'], g.loc[~worse, 'dist_ko'], s=30,
                       color=GENE_COLOR, edgecolor='white', linewidths=0.5, zorder=4,
                       label='%s KO — moved CLOSER' % lab)
        if worse.any():
            ax.scatter(g.loc[worse, 'dist_wt'], g.loc[worse, 'dist_ko'], s=30,
                       color=WORSE_COLOR, edgecolor='white', linewidths=0.5, zorder=4,
                       label='%s KO — moved AWAY' % lab)
        beat = int((g['dist_ko'] < g['dist_typical_ko']).sum())
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('distance to control, WT')
        # Three short lines, not one long one: a single-line title is wider than the panel
        # and runs into its neighbour.
        ax.set_title('%s\ndistance %+.1f%%  (typical %+.1f%%)\nbetter than typical: %d/%d'
                     % (lab, g['pct_change'].median(), g['pct_change_typical'].median(),
                        beat, len(g)), fontsize=9.5)
        ax.legend(loc='upper left', frameon=False, fontsize=6.5, markerscale=1.1)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel('distance to control, after KO')

    unit = {'gap': 'disease-gaps (|disease mean - control mean| per reaction)',
            'sd': 'control SDs', 'raw': 'raw flux units'}[scale_mode]
    devtxt = ('|v - control|' if deviation == 'abs'
              else 'max(0, disease-direction deviation)  [BIASED: credits overshoot]')
    fig.suptitle('PROTOTYPE — distance to control, before vs after KO (%s)'
                 % cohort.replace('_', '-'), fontsize=13, fontweight='medium')
    cap = ('distance(v) = sum over DE reactions of %s, in %s.\n'
           '%d reactions used; %d Status="Reversed" dropped (no defined disease direction).\n'
           'Below the dashed line = the KO moved that patient CLOSER to control. Grey = what a '
           'TYPICAL knockout does in the SAME patient; a gene must beat its own grey point.\n'
           'Unlike mTS this is NOT normalised by total change: it reports how far the patient '
           'actually ended up, not what fraction of the movement helped.\n'
           'MEASURED: every knockout INCREASES this distance, because Status gives one '
           'cohort-average direction per reaction while ~35%% of each patient\'s reactions '
           '(CD 37.7%%, UC 33.8%%) already sit on the control side —\n'
           'so the cohort direction pushes them the wrong way. mTS and recovery_proj_pca both '
           'inherit that assumption and cannot show it. Display only — [P0] open, '
           'hypothesis-only.' % (devtxt, unit, n_rxn, n_dropped))
    fig.text(0.01, -0.02, cap, fontsize=7.5, color='#666666', va='top', ha='left')

    slug = '_'.join(labels[:3]) + ('_and%dmore' % (len(labels) - 3) if len(labels) > 3 else '')
    base = os.path.join(run_dir, 'viz', 'PROTO_control_distance_%s_%s_%s%s'
                        % (cohort, slug, scale_mode,
                           '' if deviation == 'abs' else '_dis'))
    fig.savefig(base + '.png', format='png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(base + '.svg', format='svg', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return base + '.png'


def resolve_genes(viz_dir, names):
    gm = pd.read_csv(os.path.join(viz_dir, 'gene_map.csv'), dtype=str)
    id_col = gm.columns[0]
    sym2id = {}
    for _, r in gm.iterrows():
        sym2id.setdefault(str(r['symbol']).upper(), str(r[id_col]))
    known = set(gm[id_col].astype(str))
    out = []
    for nm in names:
        k = nm.strip()
        if k in known:
            out.append((k, k))
        elif k.upper() in sym2id:
            out.append((sym2id[k.upper()], k))
        else:
            raise SystemExit("FATAL gene not in this run's model: %s" % k)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--genes', default=','.join(DEFAULT_GENES))
    ap.add_argument('--scale', choices=['gap', 'sd', 'raw'], default='gap')
    ap.add_argument('--deviation', choices=['abs', 'disease'], default='abs')
    ap.add_argument('--run', action='append', metavar='LABEL=DIR')
    args = ap.parse_args()

    runs = DEFAULT_RUNS
    if args.run:
        runs = []
        for spec in args.run:
            if '=' not in spec:
                raise SystemExit('FATAL --run needs LABEL=DIR, got: %s' % spec)
            lab, _, path = spec.partition('=')
            runs.append((lab.strip(), path.strip().rstrip('/')))
    names = [g.strip() for g in args.genes.split(',') if g.strip()]

    allres = []
    for cohort, run_dir in runs:
        if not os.path.isdir(run_dir):
            raise SystemExit('FATAL not a directory: %s' % run_dir)
        print('=== %s  (%s)' % (cohort, run_dir))
        gene_ids = resolve_genes(os.path.join(run_dir, 'viz'), names)
        res, n_rxn, n_dropped, top_share = score_cohort(
            cohort, run_dir, gene_ids, args.scale, args.deviation)
        if res.empty:
            raise SystemExit('FATAL no usable patients for %s' % cohort)

        assert (res['dist_wt'] >= 0).all() and (res['dist_ko'] >= 0).all(), \
            'distance must be non-negative'
        c, s, scale, _ = disease_axis(run_dir, args.scale)
        self_d = distance_matrix(c.to_frame('c'), c, s, scale, args.deviation)[0]
        assert abs(self_d) < 1e-9, 'distance(control mean) must be 0, got %g' % self_d

        png = render(cohort, run_dir, res, args.scale, args.deviation,
                     n_rxn, n_dropped, names)
        print('    reactions used %d (dropped %d Reversed) | top reaction = %.1f%% of the '
              'WT distance' % (n_rxn, n_dropped, top_share))
        if top_share > 50:
            print('    WARNING: one reaction dominates the metric — results are about that '
                  'reaction, not the network. Try --scale gap.')
        for lab in names:
            g = res[res['gene'] == lab]
            if g.empty:
                continue
            print('    %-8s distance %+7.2f%% (typical %+6.2f%%) | closer in %2d/%d | better '
                  'than typical in %2d/%d'
                  % (lab, g['pct_change'].median(), g['pct_change_typical'].median(),
                     int((g['dist_ko'] < g['dist_wt']).sum()), len(g),
                     int((g['dist_ko'] < g['dist_typical_ko']).sum()), len(g)))
        print('    wrote %s' % png)
        allres.append(res)

    out = pd.concat(allres, ignore_index=True)
    dest = 'PROTO_control_distance_%s%s.csv' % (args.scale,
                                                '' if args.deviation == 'abs' else '_dis')
    out.to_csv(dest, index=False)
    print('\nwrote %s  (%d rows)' % (dest, len(out)))


if __name__ == '__main__':
    sys.exit(main())
