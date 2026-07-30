#!/usr/bin/env python
"""PROTOTYPE — is the [P0] "global drift" just the WT-anchor projection?

STANDALONE. Not wired into run_pipeline; modifies nothing under tf_metabolism/.
Outputs are PROTO_-prefixed so they cannot overwrite pipeline figures.

Design: docs/superpowers/specs/2026-07-31-moma-null-rebaseline-design.md

THE OBSERVATION THAT MOTIVATES THIS
-----------------------------------
Splitting every knockout into the part shared by all genes and the part specific to the
gene, for one representative patient per cohort:

                    shared by every gene | least-perturbing gene | gene-specific | ratio
    CD (H4031)                      2846 |                  2842 |          47.6 | 59.7x
    UC (H4035)                      4232 |                  4230 |          76.5 | 55.3x

Not a bookkeeping error: flux2.csv is bit-identical to MOMA's actual anchor
(flux2/Flux_prediction_<sample>.csv), max difference 0 over all 7,975 reactions. The
offset comes from the anchor itself -- the flux2 prediction leaves many reactions pinned
at +-1000 (MAR09005 -1000 -> -514.68, MAR04896 -1000 -> -577.25, ...), and MOMA's
quadratic objective shrinks them. The pull is broad: the top 5 reactions are only 11.6%
of the squared offset, and 0 of 2,172 genes reproduce the anchor.

mTS is computed on delta = KO - WT, so it CONTAINS that offset. CD's background of +0.2563
with 100% of genes positive may therefore be the projection's direction, shared by
construction rather than by biology. This script tests exactly that.

WHAT IT DOES
------------
Stage A (no solver): re-score every knockout against a NULL baseline -- the median over all
  ~2,170 genes in the same patient -- instead of against WT, using the SAME scoring code for
  both so any difference is attributable to the baseline and not to the implementation.
Stage B (Gurobi): the consensus is a proxy, since every column it averages does have some
  gene knocked out. So run the TRUE unknocked MOMA on a few patients and measure how far the
  proxy is from it, relative to the gene-specific scale.

A PASSING RESULT DOES NOT CLOSE [P0]. It would EXPLAIN the drift, which is a different and
prior thing. Still mTS-only (no bTS/wTS), still not a significance test.

Usage:
    python prototype_moma_null_rebaseline.py
    python prototype_moma_null_rebaseline.py --stage-b-patients 0     # skip the solver stage
"""
from __future__ import print_function

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from prototype_rmta_ts import desired_direction
from tf_metabolism.utils import set_publication_style


RUNS = [('CD_ileum', 'output_CD_ileum_Human2_tINIT_Human2_tINIT'),
        ('UC_rectum', 'output_UC_rectum_Human2_tINIT_Human2_tINIT')]
MARK_GENES = ['SLC7A2', 'ENO1', 'ACOX1', 'CPT2']
MIN_PATIENT_FRAC = 0.5

# Numbers this script must reproduce on the WT branch, from prototype_rmta_ts.py. If the
# regression gate fails the reimplementation is unfaithful and nothing else here means
# anything -- do NOT relax these to match the code.
EXPECTED_WT = {'CD_ileum': (+0.2563, 100.0), 'UC_rectum': (-0.0327, 1.9)}

WT_COLOR, NULL_COLOR = '#999999', '#0072B2'


def transformation_score(KO, baseline, want):
    """rMTA Eq. 8 for every column of KO against `baseline`. ONE implementation, two uses.

    TS = [sum_success |delta| - sum_unsuccess |delta|] / sum_all |delta|, bounded [-1, +1].
    `want` is +1 where the reaction must INCREASE to move disease->control.
    """
    delta = KO - baseline[:, None]
    mag = np.abs(delta)
    right = np.sign(delta) == np.sign(want)[:, None]
    den = mag.sum(axis=0)
    num = np.where(right, mag, -mag).sum(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 0, num / den, np.nan)


def stage_a(cohort, run_dir):
    """Score every knockout against BOTH baselines, per patient."""
    want = desired_direction(run_dir)
    wt_all = pd.read_csv(os.path.join(run_dir, 'flux2.csv'), index_col=0).fillna(0.0)
    files = sorted(glob.glob(os.path.join(run_dir, 'targeting_results',
                                          'MOMA_target_results_*.csv')))
    ts_wt, ts_null, offsets = {}, {}, []
    for n, f in enumerate(files, 1):
        sample = os.path.basename(f).replace('MOMA_target_results_', '').replace('.csv', '')
        if sample not in wt_all.columns:
            print('    skip %s (no WT flux column)' % sample)
            continue
        ko = pd.read_csv(f, index_col=0)
        rxns = [r for r in want.index if r in ko.index]
        if not rxns:
            print('    skip %s (no DE reactions in its MOMA output)' % sample)
            continue
        K = ko.loc[rxns].astype(np.float64).values
        w = want.reindex(rxns).values
        wt = wt_all.reindex(rxns)[sample].fillna(0.0).values
        null = np.median(K, axis=1)          # consensus over ALL genes = the proxy null

        ts_wt[sample] = pd.Series(transformation_score(K, wt, w), index=ko.columns)
        ts_null[sample] = pd.Series(transformation_score(K, null, w), index=ko.columns)

        # Offset diagnostic on the FULL reaction set, matching how it was first measured.
        rx_all = [r for r in ko.index if r in wt_all.index]
        Kf = ko.loc[rx_all].astype(np.float64).values
        wf = wt_all.reindex(rx_all)[sample].fillna(0.0).values
        Mf = np.median(Kf, axis=1)
        offsets.append(dict(cohort=cohort, sample=sample,
                            offset_consensus_wt=float(np.linalg.norm(Mf - wf)),
                            gene_specific_median=float(np.median(
                                np.linalg.norm(Kf - Mf[:, None], axis=0)))))
        print('    [%d/%d] %s: %d DE reactions | median TS  wt %+.3f -> null %+.3f'
              % (n, len(files), sample, len(rxns),
                 np.nanmedian(ts_wt[sample].values), np.nanmedian(ts_null[sample].values)))
        sys.stdout.flush()

    return pd.DataFrame(ts_wt), pd.DataFrame(ts_null), pd.DataFrame(offsets)


def summarise(TS, label):
    """Per-gene median across patients, min-patient filtered, plus background stats."""
    med = TS.median(axis=1)
    npat = TS.notna().sum(axis=1)
    keep = npat >= MIN_PATIENT_FRAC * npat.max()
    f = med[keep]
    return dict(label=label, med=med, npat=npat, keep=keep,
                background=float(f.median()), pct_pos=100.0 * float((f > 0).mean()),
                n=int(keep.sum()))


def stage_b(cohort, run_dir, n_patients):
    """True unknocked MOMA on a few patients; how far is the consensus proxy from it?"""
    if n_patients <= 0:
        return pd.DataFrame()
    try:
        from cobra.io import read_sbml_model
        from tf_metabolism.metabolic_simulation import Simulator
    except Exception as e:
        print('    Stage B skipped — cannot import solver stack: %s' % e)
        return pd.DataFrame()

    wt_all = pd.read_csv(os.path.join(run_dir, 'flux2.csv'), index_col=0).fillna(0.0)
    files = sorted(glob.glob(os.path.join(run_dir, 'targeting_results',
                                          'MOMA_target_results_*.csv')))[:n_patients]
    rows = []
    for f in files:
        sample = os.path.basename(f).replace('MOMA_target_results_', '').replace('.csv', '')
        mdl = glob.glob(os.path.join(run_dir, 'condition2', 'reconstructed_models',
                                     sample, 'functional_%s.xml' % sample))
        if not mdl:
            mdl = glob.glob(os.path.join(run_dir, 'condition2', 'reconstructed_models',
                                         '**', 'functional_%s.xml' % sample), recursive=True)
        fluxf = glob.glob(os.path.join(run_dir, 'flux2', 'Flux_prediction_%s.csv' % sample))
        if not mdl or not fluxf:
            print('    %s: model or flux file not found — skipped' % sample)
            continue

        # Rebuild MOMA's wild_flux EXACTLY as __main__.py does, including the setdefault(0.0)
        # that fills every model reaction absent from the prediction file.
        flux_dist = {}
        with open(fluxf[0]) as fp:
            for line in fp:
                p = line.strip().split(',')
                if len(p) >= 2:
                    flux_dist[p[0].strip()] = float(p[1].strip())
        cobra_model = read_sbml_model(mdl[0])
        for rxn in cobra_model.reactions:
            flux_dist.setdefault(rxn.id, 0.0)
        non_boundary = {rxn.id for rxn in cobra_model.reactions if not rxn.boundary}

        obj = Simulator.Simulator()
        obj.load_cobra_model(cobra_model)
        _, _, null_true = obj.run_MOMA(wild_flux=flux_dist, flux_constraints={})
        if null_true is False:
            print('    %s: unknocked MOMA infeasible — skipped' % sample)
            continue
        nt = pd.Series({r: v for r, v in null_true.items() if r in non_boundary})

        ko = pd.read_csv(f, index_col=0)
        rx = [r for r in ko.index if r in nt.index and r in wt_all.index]
        K = ko.loc[rx].astype(np.float64).values
        cons = np.median(K, axis=1)
        wt = wt_all.reindex(rx)[sample].fillna(0.0).values
        t = nt.loc[rx].values
        gene_specific = float(np.median(np.linalg.norm(K - cons[:, None], axis=0)))
        rows.append(dict(cohort=cohort, sample=sample,
                         consensus_vs_wt=float(np.linalg.norm(cons - wt)),
                         true_null_vs_wt=float(np.linalg.norm(t - wt)),
                         true_null_vs_consensus=float(np.linalg.norm(t - cons)),
                         gene_specific_median=gene_specific,
                         ratio_to_gene_specific=float(np.linalg.norm(t - cons)) /
                         max(gene_specific, 1e-12)))
        print('    %s: ||true_null - consensus|| %.4g  vs gene-specific %.4g  (ratio %.2f)'
              % (sample, rows[-1]['true_null_vs_consensus'], gene_specific,
                 rows[-1]['ratio_to_gene_specific']))
        sys.stdout.flush()
    return pd.DataFrame(rows)


def render(cohort, run_dir, s_wt, s_null, sym):
    set_publication_style()
    fig, ax = plt.subplots(figsize=(7.6, 5.0), facecolor='white')
    a = s_wt['med'][s_wt['keep']].values
    b = s_null['med'][s_null['keep']].values
    lo, hi = min(a.min(), b.min()), max(a.max(), b.max())
    bins = np.linspace(lo, hi, 90)
    ax.hist(a, bins=bins, color=WT_COLOR, alpha=0.75,
            label='baseline = WT   (bg %+.4f, %.0f%% positive)'
                  % (s_wt['background'], s_wt['pct_pos']))
    ax.hist(b, bins=bins, color=NULL_COLOR, alpha=0.65,
            label='baseline = MOMA null   (bg %+.4f, %.0f%% positive)'
                  % (s_null['background'], s_null['pct_pos']))
    ax.axvline(0, color='#444444', lw=1.0)
    ax.axvline(s_wt['background'], color=WT_COLOR, lw=1.4, ls='--')
    ax.axvline(s_null['background'], color=NULL_COLOR, lw=1.4, ls='--')

    ymax = ax.get_ylim()[1]
    for i, g in enumerate(MARK_GENES):
        gid = sym.get(g)
        if gid is None or gid not in s_null['med'].index:
            continue
        v = s_null['med'][gid]
        ax.annotate(g, xy=(v, ymax * 0.04), xytext=(v, ymax * (0.30 + 0.13 * (i % 3))),
                    fontsize=7.5, ha='center', color='#333333',
                    arrowprops=dict(arrowstyle='-', lw=0.5, color='#999999'))

    ax.set_xlabel('rMTA transformation score (per-gene median across patients)')
    ax.set_ylabel('genes')
    ax.set_title('PROTOTYPE — mTS background before vs after re-baselining (%s)'
                 % cohort.replace('_', '-'), fontsize=12)
    ax.legend(loc='upper left', frameon=False, fontsize=8)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    cap = ('Grey: every knockout scored against the patient WT — the baseline the pipeline and '
           'prototype_rmta_ts.py use.\n'
           'Blue: scored against the MOMA null (median over all genes in the same patient), '
           'which removes the offset every knockout shares.\n'
           'Labelled genes are positioned on the BLUE distribution. Display only — mTS is '
           'median-only (no bTS/wTS) and this does not close [P0].')
    fig.text(0.01, -0.02, cap, fontsize=7.5, color='#666666', va='top', ha='left')
    base = os.path.join(run_dir, 'viz', 'PROTO_moma_null_rebaseline_%s' % cohort)
    fig.savefig(base + '.png', format='png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(base + '.svg', format='svg', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return base + '.png'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--stage-b-patients', type=int, default=5)
    args = ap.parse_args()

    per_gene, validation, failures = [], [], []
    for cohort, run_dir in RUNS:
        print('=== %s  (%s)' % (cohort, run_dir))
        TSW, TSN = stage_a(cohort, run_dir)[:2]
        s_wt, s_null = summarise(TSW, 'wt'), summarise(TSN, 'null')

        # ---- Regression gate -------------------------------------------------
        exp_bg, exp_pos = EXPECTED_WT[cohort]
        ok = abs(s_wt['background'] - exp_bg) < 5e-4 and abs(s_wt['pct_pos'] - exp_pos) < 0.6
        print('    REGRESSION (WT branch): background %+.4f (expected %+.4f) | %.1f%% positive '
              '(expected %.1f%%) -> %s'
              % (s_wt['background'], exp_bg, s_wt['pct_pos'], exp_pos, 'PASS' if ok else 'FAIL'))
        if not ok:
            failures.append(cohort)

        print('    RE-BASELINED:           background %+.4f | %.1f%% positive  (n=%d genes)'
              % (s_null['background'], s_null['pct_pos'], s_null['n']))

        both = s_wt['keep'] & s_null['keep']
        rho, p = spearmanr(s_wt['med'][both], s_null['med'][both])
        print('    ranking preserved? Spearman(mTS_wt, mTS_null) = %+.3f (p=%.2g, n=%d)'
              % (rho, p, int(both.sum())))

        gm = pd.read_csv(os.path.join(run_dir, 'viz', 'gene_map.csv'), dtype=str)
        sym2id = dict(zip(gm['symbol'], gm.iloc[:, 0]))
        id2sym = dict(zip(gm.iloc[:, 0], gm['symbol']))
        ref_w = np.sort(s_wt['med'][s_wt['keep']].values)
        ref_n = np.sort(s_null['med'][s_null['keep']].values)
        print('    %-8s %>10s %>10s %>9s %>9s' % ('gene', 'mTS_wt', 'mTS_null', 'pct_wt', 'pct_null')
              if False else '    gene       mTS_wt    mTS_null    pct_wt   pct_null')
        for g in MARK_GENES:
            gid = sym2id.get(g)
            if gid is None or gid not in s_null['med'].index:
                continue
            vw, vn = s_wt['med'][gid], s_null['med'][gid]
            print('    %-8s %+9.4f  %+9.4f   %7.1f   %8.1f'
                  % (g, vw, vn,
                     100 * np.searchsorted(ref_w, vw) / len(ref_w),
                     100 * np.searchsorted(ref_n, vn) / len(ref_n)))

        png = render(cohort, run_dir, s_wt, s_null, sym2id)
        print('    wrote %s' % png)

        t = pd.DataFrame({'Gene': s_null['med'].index,
                          'Symbol': [id2sym.get(g, g) for g in s_null['med'].index],
                          'mTS_wt': s_wt['med'].values, 'mTS_null': s_null['med'].values,
                          'n_patients': s_null['npat'].values})
        t.insert(0, 'cohort', cohort)
        per_gene.append(t)

        print('    --- Stage B: true unknocked MOMA on %d patients' % args.stage_b_patients)
        v = stage_b(cohort, run_dir, args.stage_b_patients)
        if not v.empty:
            validation.append(v)
            print('    proxy error / gene-specific scale: median %.2f'
                  % v['ratio_to_gene_specific'].median())
        print()

    pd.concat(per_gene, ignore_index=True).to_csv('PROTO_moma_null_rebaseline.csv', index=False)
    print('wrote PROTO_moma_null_rebaseline.csv')
    if validation:
        pd.concat(validation, ignore_index=True).to_csv('PROTO_moma_null_validation.csv',
                                                        index=False)
        print('wrote PROTO_moma_null_validation.csv')
    if failures:
        print('\nREGRESSION GATE FAILED for: %s — results above are NOT interpretable.'
              % ', '.join(failures))
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
