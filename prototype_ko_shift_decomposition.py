#!/usr/bin/env python
"""PROTOTYPE — is a gene's UMAP knockout shift SPECIFIC, or just the cohort-wide drift?

STANDALONE. Not wired into run_pipeline; modifies nothing under tf_metabolism/.
Outputs are PROTO_-prefixed so they cannot overwrite pipeline figures.

THE QUESTION
------------
"If I knock out SLC7A2 in CD, will I see a shift in the UMAP?"

Yes — a large one. That is exactly the problem. In CD-ileum ~79% of ALL genes move
toward control and 100% score positive on the rMTA mTS, with a background median of
+0.2563. SLC7A2's mTS is +0.2878, so ~89% of its score is the cohort-wide baseline that
EVERY knockout gets. A big arrow toward the control cloud is the default behaviour of
this cohort, not evidence about this gene.

So "does it shift?" is the wrong question. The right one is:

    Does this knockout shift DIFFERENTLY from a typical knockout in the same patient?

THE DECOMPOSITION
-----------------
For one patient, every gene g has a UMAP displacement d_g = KO_g - WT (the WT point is
shared by all genes in that patient, so all displacements start from the same origin).

    D       = componentwise MEDIAN of d_g over all genes   -> the cohort-wide drift
    dev_g   = || d_g - D ||                                -> departure from typical
    par_g   = d_g . D_hat                                  -> movement along the drift

Panel A reports ||D|| and the share of a gene's shift lying along it — i.e. how much of
the visible arrow is cohort default rather than the gene.

WHAT GEOMETRY CANNOT DO (measured, not assumed)
-----------------------------------------------
The first version of panel B scored genes by ||d_g - D||, the departure from the drift,
as a per-patient percentile. Tested over all ~2,200 genes it carries essentially NO
information about whether the knockout helps:

    Spearman(departure-from-drift, mTS excess) = +0.006 (CD, p=0.78) / +0.159 (UC)
    top-5% vs bottom-5% genes by mTS excess:  4.5 (CD) / 0.1 (UC) percentile points apart

UC's ACOX1 (mTS 99.95th pctile) and CPT2 (7.9th) sat on top of each other at 90.9 vs 86.3.
Even the SIGNED recovery_proj_umap manages only +0.087 / +0.203.

The reason is structural: mTS is per-reaction and direction-aware in flux space, while the
UMAP is a 2D nonlinear embedding carrying no per-reaction sign. So panel B now reads the
flux-space quantity directly — the rMTA TS ranked against all genes in the SAME patient.
Departure-from-drift is retained in the CSV as a specificity measure, but it is NOT a
discriminator and is labelled as such wherever it is printed.

Neither panel is a significance test; [P0] stays open and mTS is median-only (no bTS/wTS).

Nothing is re-simulated. Reads the stored per-patient UMAP MOMA projections, so it needs
a run whose projection came from the current predict-frame moma_project.R (the two Human2
runs; the Recon2M.2 runs lack wt_x/wt_y and are out of scope).

Usage:
    python prototype_ko_shift_decomposition.py                    # CD: SLC7A2 vs ENO1
    python prototype_ko_shift_decomposition.py --run UC=output_... --genes ACOX1,CPT2
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


DEFAULT_RUN = ('CD_ileum', 'output_CD_ileum_Human2_tINIT_Human2_tINIT')
DEFAULT_GENES = ['SLC7A2', 'ENO1']

GENE_COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7']
DRIFT_COLOR = '#333333'
CLOUD_COLOR = '#CCCCCC'


def load_ts_percentiles(viz_dir, gene_ids):
    """Per-patient rMTA TS percentile of each gene among ALL genes in that same patient.

    THIS is the discriminating quantity, and it is deliberately NOT geometric. An earlier
    version of panel B scored genes by how far their UMAP displacement departed from the
    cohort drift; tested over all ~2,200 genes that turned out to carry essentially no
    information about mTS -- Spearman +0.006 (CD) / +0.159 (UC), and the top-5% and
    bottom-5% genes by mTS excess were separated by 4.5 and 0.1 percentile points
    respectively. UC's ACOX1 (mTS 99.95th pctile) and CPT2 (7.9th) landed on top of each
    other at 90.9 vs 86.3.

    The reason is structural: mTS is per-reaction and direction-aware in flux space, while
    the UMAP is a 2D nonlinear embedding carrying no per-reaction sign. No geometric
    summary of the embedding can recover it -- even the signed recovery_proj_umap manages
    only +0.087 / +0.203. So the panel has to read the flux-space quantity directly.
    """
    path = os.path.join(viz_dir, 'prototype_rmta_TS_matrix.csv')
    if not os.path.exists(path):
        raise SystemExit(
            'FATAL no per-patient TS matrix: %s\n'
            '      Generate it first:  python prototype_rmta_ts.py %s\n'
            '      (~10 min; re-reads the stored MOMA matrices.)'
            % (path, os.path.dirname(viz_dir.rstrip('/'))))

    TS = pd.read_csv(path, index_col=0)
    # Rank within each patient (column): "does this KO beat a typical KO in THIS patient".
    pct = TS.rank(axis=0, pct=True) * 100.0
    rows = []
    for gid, label in gene_ids:
        if gid not in TS.index:
            continue
        for sample in TS.columns:
            t, p = TS.loc[gid, sample], pct.loc[gid, sample]
            if np.isfinite(t) and np.isfinite(p):
                rows.append(dict(sample=str(sample), gene=label,
                                 ts=float(t), ts_pctile=float(p)))
    return pd.DataFrame(rows)


def gene_slug(labels):
    """Filename fragment naming the plotted genes.

    Every output is named after the genes it actually contains. Without this, the CD and
    UC figures differ only by a cohort suffix while plotting entirely different genes —
    which is exactly how a UC (ACOX1/CPT2) figure gets mistaken for the CD (SLC7A2/ENO1)
    one. Capped so a long --genes list cannot produce an unwieldy path.
    """
    if not labels:
        return 'none'
    if len(labels) > 3:
        return '%s_and%dmore' % (labels[0], len(labels) - 1)
    return '_'.join(labels)


def resolve_genes(viz_dir, names):
    """Map HGNC symbols (or raw model ids) to model gene ids. Fails loud on a miss."""
    gm_path = os.path.join(viz_dir, 'gene_map.csv')
    if not os.path.exists(gm_path):
        raise SystemExit('FATAL no gene map: %s' % gm_path)
    gm = pd.read_csv(gm_path)
    # The header says "entrez" but the column holds ENSG ids for Human2 runs (CLAUDE.md).
    id_col = gm.columns[0]
    sym2id = {}
    for _, r in gm.iterrows():
        sym2id.setdefault(str(r['symbol']).upper(), str(r[id_col]))
    known_ids = set(gm[id_col].astype(str))

    out = []
    for nm in names:
        key = nm.strip()
        if key in known_ids:
            out.append((key, key))
        elif key.upper() in sym2id:
            out.append((sym2id[key.upper()], key))
        else:
            raise SystemExit('FATAL gene not in this run\'s model: %s' % key)
    return out


def per_patient_decomposition(viz_dir, gene_ids):
    """Per patient: drift vector, and each gene's departure from it + its percentile."""
    files = sorted(glob.glob(os.path.join(
        viz_dir, 'UMAP_Results_MOMA_target_results_*.csv')))
    if not files:
        raise SystemExit('FATAL no UMAP_Results_MOMA_target_results_*.csv in %s' % viz_dir)

    recs = []
    per_patient_raw = {}
    for f in files:
        sample = os.path.basename(f).replace(
            'UMAP_Results_MOMA_target_results_', '').replace('.csv', '')
        d = pd.read_csv(f, index_col=0)
        need = ('X', 'Y', 'wt_x', 'wt_y')
        if any(c not in d.columns for c in need):
            raise SystemExit(
                'FATAL %s lacks %s — this run predates the predict-frame moma_project.R '
                '(commit e388856). Only the Human2 runs are supported.'
                % (f, ', '.join(need)))

        dx = d['X'].values - d['wt_x'].values
        dy = d['Y'].values - d['wt_y'].values
        # Drift = what a knockout does by default in this patient. Median, not mean:
        # a handful of extreme knockouts would otherwise define "typical".
        Dx, Dy = np.median(dx), np.median(dy)
        dev = np.sqrt((dx - Dx) ** 2 + (dy - Dy) ** 2)      # departure from typical
        mag = np.sqrt(dx ** 2 + dy ** 2)                     # total displacement
        Dn = np.sqrt(Dx ** 2 + Dy ** 2)
        par = (dx * Dx + dy * Dy) / Dn if Dn > 0 else np.zeros_like(dx)

        proj = (d['recovery_proj_umap'].values if 'recovery_proj_umap' in d.columns
                else np.full(len(d), np.nan))
        dev_sorted = np.sort(dev)
        proj_sorted = np.sort(proj[np.isfinite(proj)])

        per_patient_raw[sample] = dict(
            d=d, dx=dx, dy=dy, D=(Dx, Dy), dev=dev, drift_norm=Dn,
            wt=(float(d['wt_x'].iloc[0]), float(d['wt_y'].iloc[0])))

        for gid, label in gene_ids:
            if gid not in d.index:
                continue
            i = d.index.get_loc(gid)
            if isinstance(i, slice) or isinstance(i, np.ndarray):
                i = np.arange(len(d))[i][0]
            recs.append(dict(
                sample=sample, gene=label, gene_id=gid,
                drift_norm=Dn,
                shift_norm=float(mag[i]),
                drift_share=float(par[i]) / float(mag[i]) if mag[i] > 0 else np.nan,
                dev=float(dev[i]),
                dev_pctile=100.0 * np.searchsorted(dev_sorted, dev[i]) / len(dev_sorted),
                recovery_proj_umap=float(proj[i]),
                proj_pctile=(100.0 * np.searchsorted(proj_sorted, proj[i]) /
                             len(proj_sorted) if len(proj_sorted) else np.nan),
                toward=bool(proj[i] > 0) if np.isfinite(proj[i]) else False,
            ))

    return pd.DataFrame(recs), per_patient_raw


def render(cohort, viz_dir, res, raw, gene_ids):
    """Panel A: one patient's KO cloud with drift + gene arrows. Panel B: dev percentiles."""
    # Representative patient = median drift magnitude, so the panel is neither the most
    # nor the least drifted patient.
    order = sorted(raw.items(), key=lambda kv: kv[1]['drift_norm'])
    sample, r = order[len(order) // 2]

    set_publication_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.0, 5.8), facecolor='white',
                                   gridspec_kw={'width_ratios': [1.25, 1.0]})

    # ---- Panel A -------------------------------------------------------------
    ref_path = os.path.join(viz_dir, 'UMAP_predict_ref_df.csv')
    if os.path.exists(ref_path):
        ref = pd.read_csv(ref_path, index_col=0)
        ctrl = ref[ref['condition'].astype(str).str.lower().str.contains('control')]
        axA.scatter(ctrl['UMAP1'], ctrl['UMAP2'], s=26, color='#9ECAE1',
                    edgecolor='white', linewidths=0.5, zorder=2, label='nonIBD control')

    d = r['d']
    axA.scatter(d['X'], d['Y'], s=5, color=CLOUD_COLOR, linewidths=0, zorder=1,
                label='every gene KO (n=%d)' % len(d))
    wx, wy = r['wt']
    axA.scatter([wx], [wy], s=90, marker='*', color='#000000', zorder=6,
                label='patient WT')
    Dx, Dy = r['D']
    axA.annotate('', xy=(wx + Dx, wy + Dy), xytext=(wx, wy), zorder=5,
                 arrowprops=dict(arrowstyle='-|>', lw=2.2, color=DRIFT_COLOR))
    axA.scatter([wx + Dx], [wy + Dy], s=45, marker='X', color=DRIFT_COLOR, zorder=6,
                label='cohort drift (median KO)')

    for (gid, label), col in zip(gene_ids, GENE_COLORS):
        if gid not in d.index:
            continue
        gx, gy = float(d.loc[gid, 'X']), float(d.loc[gid, 'Y'])
        axA.annotate('', xy=(gx, gy), xytext=(wx, wy), zorder=7,
                     arrowprops=dict(arrowstyle='-|>', lw=1.6, color=col))
        axA.scatter([gx], [gy], s=55, color=col, edgecolor='white', linewidths=0.7,
                    zorder=8, label='%s KO' % label)

    axA.set_title('A — patient %s (median-drift patient)' % sample, fontsize=11)
    axA.set_xlabel('UMAP1')
    axA.set_ylabel('UMAP2')
    axA.legend(loc='best', frameon=False, fontsize=7.5, markerscale=1.1)
    for s in ('top', 'right'):
        axA.spines[s].set_visible(False)

    # ---- Panel B: does the KO beat a typical KO at moving THIS patient toward control? --
    labels = [lab for _, lab in gene_ids if lab in set(res['gene'])]
    axB.axhline(50, color='#888888', lw=1.0, ls='--', zorder=1)
    axB.text(0.985, 51, 'a typical knockout', transform=axB.get_yaxis_transform(),
             ha='right', va='bottom', fontsize=7.5, color='#666666')
    have_ts = 'ts_pctile' in res.columns and res['ts_pctile'].notna().any()
    for i, lab in enumerate(labels):
        v = res[res['gene'] == lab]['ts_pctile'].dropna().values
        if not len(v):
            continue
        x = np.full(len(v), i, dtype=float) + np.linspace(-0.16, 0.16, len(v))
        axB.scatter(x, v, s=26, color=GENE_COLORS[i % len(GENE_COLORS)],
                    edgecolor='white', linewidths=0.5, zorder=3)
        axB.plot([i - 0.30, i + 0.30], [np.median(v)] * 2, color='#000000', lw=2.0,
                 zorder=4)
        axB.annotate('median %.0f' % np.median(v), xy=(i + 0.33, np.median(v)),
                     fontsize=8, va='center', ha='left', color='#333333')
    axB.set_xticks(range(len(labels)))
    axB.set_xticklabels(labels)
    axB.set_ylim(0, 105)
    axB.set_ylabel('rMTA transformation-score percentile\n'
                   '(among all genes, within each patient)')
    axB.set_title('B — does the KO beat a typical KO? one point = one patient',
                  fontsize=11)
    if not have_ts:
        axB.text(0.5, 0.5, 'no per-patient TS matrix', transform=axB.transAxes,
                 ha='center', va='center', color='#999999')
    for s in ('top', 'right'):
        axB.spines[s].set_visible(False)

    # The genes MUST be in the title: two cohorts' figures are otherwise identical apart
    # from a directory, and the gene names appear only in a legend entry and a tick label.
    fig.suptitle('PROTOTYPE — KO shift decomposition: %s  (%s)'
                 % (' vs '.join(labels), cohort.replace('_', '-')),
                 fontsize=13, fontweight='medium')
    cap = ('A: every gene KO in a patient starts from the same WT point; the cohort drift is the '
           'MEDIAN displacement over all genes — what a knockout does by default.\n'
           'B: the rMTA transformation score (per-reaction, direction-aware), ranked against all '
           'genes in the SAME patient. ~50 = no better than a typical KO.\n'
           'The two panels answer different questions on purpose: A is how much of the visible '
           'movement is cohort-wide drift, B is whether the movement helps.\n'
           'Geometry alone cannot answer B — departure-from-drift correlates with mTS at only '
           '+0.006 (CD) / +0.159 (UC), so it was replaced by this panel.\n'
           'Display only — [P0] competitive null still open, hypothesis-only. mTS is median-only, '
           'no bTS/wTS.')
    fig.text(0.01, -0.02, cap, fontsize=7.5, color='#666666', va='top', ha='left')

    base = os.path.join(viz_dir, 'PROTO_ko_shift_decomposition_%s_%s'
                                 % (cohort, gene_slug(labels)))
    fig.savefig(base + '.png', format='png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(base + '.svg', format='svg', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return base + '.png'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', metavar='LABEL=DIR',
                    default='%s=%s' % DEFAULT_RUN)
    ap.add_argument('--genes', default=','.join(DEFAULT_GENES))
    ap.add_argument('--out-csv', default=None)
    args = ap.parse_args()

    if '=' not in args.run:
        raise SystemExit('FATAL --run needs LABEL=DIR, got: %s' % args.run)
    cohort, _, run_dir = args.run.partition('=')
    cohort, run_dir = cohort.strip(), run_dir.strip().rstrip('/')
    viz_dir = os.path.join(run_dir, 'viz')
    if not os.path.isdir(viz_dir):
        raise SystemExit('FATAL not a directory: %s' % viz_dir)

    gene_ids = resolve_genes(viz_dir, args.genes.split(','))
    res, raw = per_patient_decomposition(viz_dir, gene_ids)
    if res.empty:
        raise SystemExit('FATAL none of the requested genes appear in any patient model.')

    ts = load_ts_percentiles(viz_dir, gene_ids)
    res = res.merge(ts, on=['sample', 'gene'], how='left')

    drift = pd.Series({s: r['drift_norm'] for s, r in raw.items()})
    print('=== %s  (%d patients)' % (cohort, len(raw)))
    print('    cohort drift ||D||        : median %.3f  (range %.3f-%.3f) UMAP units'
          % (drift.median(), drift.min(), drift.max()))
    print()
    for _, lab in gene_ids:
        g = res[res['gene'] == lab]
        if g.empty:
            print('    %-8s : absent from every patient model' % lab)
            continue
        print('    %-8s present in %d/%d patients' % (lab, len(g), len(raw)))
        print('        total shift ||d||      : median %.3f  (drift is %.3f)'
              % (g['shift_norm'].median(), drift.median()))
        print('        share along the drift  : median %.0f%%'
              % (100 * g['drift_share'].median()))
        print('        moved TOWARD control   : %d/%d patients'
              % (int(g['toward'].sum()), len(g)))
        # The discriminator. Geometry (dev_pctile below) does NOT track this -- see
        # load_ts_percentiles.
        tsv = g['ts_pctile'].dropna()
        if len(tsv):
            print('    ->  rMTA TS percentile within patient : median %.1f  (%d/%d patients '
                  'above 50)' % (tsv.median(), int((tsv > 50).sum()), len(tsv)))
        else:
            print('    ->  rMTA TS percentile within patient : unavailable')
        print('        [geometry, non-discriminating] departure from drift %.3f, pctile %.1f'
              % (g['dev'].median(), g['dev_pctile'].median()))
        print()

    png = render(cohort, viz_dir, res, raw, gene_ids)
    print('    wrote %s' % png)
    present = [lab for _, lab in gene_ids if lab in set(res['gene'])]
    out_csv = args.out_csv or ('PROTO_ko_shift_decomposition_%s_%s.csv'
                               % (cohort, gene_slug(present)))
    res.to_csv(out_csv, index=False)
    print('    wrote %s  (%d rows)' % (out_csv, len(res)))


if __name__ == '__main__':
    sys.exit(main())
