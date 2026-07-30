#!/usr/bin/env python
"""PROTOTYPE — rMTA-style transformation score (mTS) on existing MOMA results.

Implements the Transformation Score of Valcarcel et al. 2019 (Bioinformatics,
"rMTA: robust metabolic transformation analysis"), Eq. 8:

    TS = [ sum_{i in Rsuccess}   |v_ref,i - v_res,i|
         - sum_{i in Runsuccess} |v_ref,i - v_res,i| ]  /  sum_{i in Rs} |v_ref,i - v_res,i|

  v_ref  = the source state          -> this patient's wild-type flux (flux2.csv)
  v_res  = the result state          -> that gene's MOMA knockout flux
  Rs     = reactions that SHOULD change -> the DE-flux reactions
  success/unsuccess = whether the knockout moved the reaction in the direction
           that takes disease toward control.

Desired direction comes from Differential_fluxes.csv `Status`: a reaction UP in
disease must come DOWN, and vice versa.

WHY THIS IS DIFFERENT from the pipeline's current recovery_proj_pca:
  * BOUNDED in [-1, +1] and normalised by the total change, so it cannot be won
    by sheer perturbation magnitude (the open [P0] problem).
  * DIRECTION-SPECIFIC per reaction: moving the wrong way is subtracted, rather
    than being invisible because it was orthogonal to a single PCA axis.

WHAT THIS IS NOT:
  * This is only rMTA's mTS (the MOMA-based score). The full rMTA criterion also
    needs bTS (best case) and wTS (worst case), which require new MTA MILP runs
    per gene per patient. Not computable from stored MOMA output.
  * Not a significance test. The percentile column is a competitive comparison
    against all other genes in the same patient (CMap-tau style), which is what
    tells you whether a gene beats a typical knockout.

Usage:
    python prototype_rmta_ts.py <run_folder> [--max-patients N] [--out FILE]
"""
from __future__ import print_function

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd


GENES_OF_INTEREST = ['ELOVL5', 'ELOVL1', 'TECR', 'ELOVL4', 'HIGD2A',
                     'HSD17B12', 'ELOVL7', 'ELOVL6', 'ELOVL2']


def desired_direction(run):
    """Per DE reaction: +1 if it must INCREASE to move disease->control, -1 if decrease.

    Differential_fluxes.csv `Status` is the direction in condition2 (disease)
    relative to condition1 (control). To move toward control we must reverse it.
    """
    df = pd.read_csv(os.path.join(run, 'Differential_fluxes.csv'), index_col=0)
    status = df['Status'].astype(str).str.upper()
    want = pd.Series(np.where(status == 'UP', -1.0, 1.0), index=df.index)
    return want[status.isin(['UP', 'DOWN'])]


def transformation_scores(run, max_patients=None):
    want = desired_direction(run)
    wt_all = pd.read_csv(os.path.join(run, 'flux2.csv'), index_col=0).fillna(0.0)

    files = sorted(glob.glob(os.path.join(run, 'targeting_results',
                                          'MOMA_target_results_*.csv')))
    if max_patients:
        files = files[:int(max_patients)]

    per_patient = {}
    for n, f in enumerate(files, 1):
        sample = os.path.basename(f).replace('MOMA_target_results_', '').replace('.csv', '')
        if sample not in wt_all.columns:
            print('  skip %s (no WT flux column)' % sample)
            continue

        # reactions x genes, restricted to the DE reactions present in this model
        ko = pd.read_csv(f, index_col=0)
        rxns = [r for r in want.index if r in ko.index]
        if not rxns:
            print('  skip %s (no DE reactions in its MOMA output)' % sample)
            continue

        ko = ko.loc[rxns].astype(np.float64)
        wt = wt_all.reindex(rxns)[sample].fillna(0.0).values[:, None]
        d = want.reindex(rxns).values[:, None]          # desired direction per reaction

        delta = ko.values - wt                          # change caused by each knockout
        mag = np.abs(delta)
        moved_right = (np.sign(delta) == np.sign(d))

        denom = mag.sum(axis=0)
        num = np.where(moved_right, mag, -mag).sum(axis=0)
        with np.errstate(invalid='ignore', divide='ignore'):
            ts = np.where(denom > 0, num / denom, np.nan)

        per_patient[sample] = pd.Series(ts, index=ko.columns)
        print('  [%d/%d] %s: %d DE reactions, %d genes, median TS %+.3f'
              % (n, len(files), sample, len(rxns), ko.shape[1],
                 np.nanmedian(ts)))
        sys.stdout.flush()

    return pd.DataFrame(per_patient)          # genes x patients


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run')
    ap.add_argument('--max-patients', type=int, default=None)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()

    print('=== rMTA-style mTS prototype: %s ===' % a.run)
    TS = transformation_scores(a.run, a.max_patients)
    if TS.empty:
        print('No usable patients.')
        return

    # Competitive background: where does each gene sit among ALL genes of the
    # SAME patient? This is the comparison the current pipeline never makes.
    pct = TS.rank(axis=0, pct=True) * 100.0

    # Persist the full genes x patients matrix, not just its per-gene medians. Collapsing
    # to a median discards the per-patient competitive comparison that `pct` computes
    # here, and recovering it costs a full ~10 min re-read of the stored MOMA matrices.
    # Downstream figures that ask "does this KO beat a typical KO IN THIS PATIENT" need it.
    mat_dest = os.path.join(a.run, 'viz', 'prototype_rmta_TS_matrix.csv')
    TS.to_csv(mat_dest)
    print('Wrote %s  (%d genes x %d patients)' % (mat_dest, TS.shape[0], TS.shape[1]))

    gm = pd.read_csv(os.path.join(a.run, 'viz', 'gene_map.csv'), dtype=str)
    id2sym = dict(zip(gm.iloc[:, 0], gm.iloc[:, 1]))

    out = pd.DataFrame({
        'Gene': TS.index,
        'Symbol': [id2sym.get(g, g) for g in TS.index],
        'mTS_median': TS.median(axis=1).values,
        'mTS_mean': TS.mean(axis=1).values,
        'pct_of_genes_beaten_median': pct.median(axis=1).values,
        'n_patients_TS_positive': (TS > 0).sum(axis=1).values,
        'n_patients': TS.notna().sum(axis=1).values,
    }).sort_values('mTS_median', ascending=False)

    dest = a.out or os.path.join(a.run, 'viz', 'prototype_rmta_mTS.csv')
    out.to_csv(dest, index=False)
    print('\nWrote %s' % dest)

    print('\n--- BACKGROUND over all %d genes ---' % len(out))
    print('   median mTS across genes: %+.3f' % out['mTS_median'].median())
    print('   genes with median mTS > 0: %d/%d (%.0f%%)'
          % ((out['mTS_median'] > 0).sum(), len(out),
             100.0 * (out['mTS_median'] > 0).mean()))

    print('\n--- TOP 10 genes by median mTS ---')
    for _, r in out.head(10).iterrows():
        print('   %-10s mTS %+.3f | beats %.0f%% of genes | TS>0 in %d/%d patients'
              % (r['Symbol'], r['mTS_median'], r['pct_of_genes_beaten_median'],
                 r['n_patients_TS_positive'], r['n_patients']))

    print('\n--- GENES OF INTEREST ---')
    sel = out[out['Symbol'].isin(GENES_OF_INTEREST)]
    for _, r in sel.iterrows():
        print('   %-10s mTS %+.3f | beats %.0f%% of genes | TS>0 in %d/%d patients'
              % (r['Symbol'], r['mTS_median'], r['pct_of_genes_beaten_median'],
                 r['n_patients_TS_positive'], r['n_patients']))


if __name__ == '__main__':
    main()
