import json
import tempfile
import argparse
import logging
import os
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import umap
import joblib
import glob

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from cobra.io import read_sbml_model, write_sbml_model

# ---------------------------------------------------------------------------
# R visualization scripts live under tf_metabolism/r_based_viz/
# ---------------------------------------------------------------------------
R_VIZ_DIR = os.path.join(os.path.dirname(__file__), 'r_based_viz')


def _r_script(name):
    """Absolute path to a bundled R visualization script."""
    return os.path.join(R_VIZ_DIR, name)


# ---------------------------------------------------------------------------
# Rscript resolver
# ---------------------------------------------------------------------------
def _find_rscript():
    """Return the Rscript executable path.

    Search order:
    1. RSCRIPT env variable (user override)
    2. shutil.which('Rscript') — searches PATH
    3. Raises RuntimeError if not found
    """
    env = os.environ.get('RSCRIPT')
    if env:
        return env
    found = shutil.which('Rscript')
    if found:
        return found
    raise RuntimeError(
        'Rscript not found. Set the RSCRIPT environment variable or add R to PATH.'
    )


def _effective_gene_map(output_viz_dir=None):
    """Path to the gene-id -> HGNC symbol map (columns: entrez,symbol). Prefers a
    run-local ``<viz>/gene_map.csv`` (written for Human-GEM: ENSG->symbol in the same
    format), else the bundled Recon2M.2 Entrez map. So both the Python readers and
    moma_viz.R resolve symbols for whichever model produced the run, unchanged."""
    if output_viz_dir:
        local = os.path.join(output_viz_dir, 'gene_map.csv')
        if os.path.exists(local):
            return local
    return os.path.join(os.path.dirname(__file__), 'data', 'metabolic_genes.csv')


def run_pathway_enrichment_plot(output_dir, output_viz_dir,
                                fdr_cutoff=0.1, title=None, r_script=None):
    """Visualize pathway enrichment results via pathway_enrichment_plot.R.

    Reads Up_regulated_pathways.csv / Down_regulated_pathways.csv from
    ``output_dir`` and writes ``{output_viz_dir}/pathway_barplot.png`` (+ .svg).

    Parameters
    ----------
    output_dir     : directory containing the enrichment CSVs
    output_viz_dir : where the figures are written
    fdr_cutoff     : Adjusted P-value threshold (default 0.1)
    title          : optional plot title suffix (e.g. disease label)
    r_script       : path to pathway_enrichment_plot.R (auto-detected if None)
    """
    import subprocess

    if r_script is None:
        r_script = _r_script('pathway_enrichment_plot.R')
    if not os.path.exists(r_script):
        raise FileNotFoundError(f'R script not found: {r_script}')

    os.makedirs(output_viz_dir, exist_ok=True)

    cmd = [
        _find_rscript(), r_script,
        '--output-dir',     output_dir,
        '--output-viz-dir', output_viz_dir,
        '--fdr-cutoff',     str(fdr_cutoff),
    ]
    if title:
        cmd += ['--title', title]

    logging.info('[R] Running: %s', ' '.join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError('pathway_enrichment_plot.R exited with code %d' % result.returncode)
    logging.info('[R] Pathway enrichment plot done.')


def visualize_pathway_enrichment(output_dir, output_viz_dir, use_r=True,
                                 fdr_cutoff=0.1, title=None):
    """Bar plot of enriched metabolic pathways, with a selectable backend.

    Mirrors the UMAP visualization style: pass ``use_r=True`` to render with the
    R script (``pathway_enrichment_plot.R``, ggplot2) or ``use_r=False`` to render
    with matplotlib. Both backends read Up/Down_regulated_pathways.csv from
    ``output_dir`` and write ``{output_viz_dir}/pathway_barplot.png`` (+ .svg).

    Parameters
    ----------
    output_dir     : directory containing the enrichment CSVs
    output_viz_dir : where the figures are written
    use_r          : True → R/ggplot2 backend, False → matplotlib backend
    fdr_cutoff     : Adjusted P-value threshold (default 0.1)
    title          : optional plot title suffix (e.g. disease label)
    """
    if use_r:
        run_pathway_enrichment_plot(output_dir, output_viz_dir,
                                    fdr_cutoff=fdr_cutoff, title=title)
    else:
        _pathway_barplot_matplotlib(output_dir, output_viz_dir,
                                    fdr_cutoff=fdr_cutoff, title=title)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
CONDITION_COLORS = {'Control': 'royalblue', 'Test': '#DC143C'}
DATASET_MARKERS  = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', 'X', '8']


def argument_parser(version=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory")
    parser.add_argument('-c1', '--input_omics_file1', required=True,
                        help="Input omics file for condition 1")
    parser.add_argument('-c2', '--input_omics_file2', required=True,
                        help="Input omics file for condition 2")
    parser.add_argument('--flux_pval', type=float, default=0.05,
                        help="DE-flux p-value cutoff")
    parser.add_argument('--flux_filter', choices=['raw_p', 'FDR', 'Bonferroni', 'None'],
                        default='raw_p', help="DE-flux p-value filter mode")
    parser.add_argument('--flux_rbc', type=float, default=0.5,
                        help="DE-flux effect-size cutoff |rank-biserial|")
    parser.add_argument('--target_fdr', type=float, default=0.05,
                        help="Cohort target-ranking BH-FDR cutoff")
    parser.add_argument('--target_effect_min', type=float, default=0.1,
                        help="Cohort target-ranking minimum effect size (recurrence proportion)")
    parser.add_argument('--target_topk', type=int, default=10,
                        help="Cohort target-ranking top-k (recurrence / bootstrap selection)")
    parser.add_argument('--target_method', choices=['auto', 'rra', 'mixedlm'], default='auto',
                        help="Cohort target-ranking inference: 'auto' (size-adaptive: RRA + "
                             "permutation always, add mixed model at n>=target_switch_n), "
                             "'rra' (RRA only), or 'mixedlm' (force mixed model)")
    parser.add_argument('--target_switch_n', type=int, default=40,
                        help="Patient-count threshold at which 'auto' adds the mixed-effects model")
    parser.add_argument('--target_nperm', type=int, default=1000,
                        help="Cohort target-ranking permutation-null resamples")
    parser.add_argument('--target_nboot', type=int, default=1000,
                        help="Cohort target-ranking bootstrap resamples")
    parser.add_argument('--model', nargs='+', default=['Recon2M.2'],
                        help="Base GEM(s), one or more (space- or comma-separated). "
                             "Choices: Recon2M.2, Human1 (Human-GEM v1.0.0), Human2 "
                             "(Human-GEM v2.0.0). Default: Recon2M.2")
    parser.add_argument('--mem_method', nargs='+', default=['tINIT'],
                        help="Model-extraction method(s), one or more (space- or "
                             "comma-separated). Choices: MBA, GIMME, FASTCORE, iMAT, "
                             "INIT, tINIT, mCADRE. Default: tINIT (only tINIT is "
                             "runnable in this environment; others fail loud).")
    return parser


def update_cobra_model(cobra_model):
    temp = tempfile.NamedTemporaryFile(prefix='temp_cobra_', suffix='.xml', dir='./', delete=True)
    temp_outfile = temp.name
    temp.close()

    write_sbml_model(cobra_model, temp_outfile)
    cobra_model = read_sbml_model(temp_outfile)
    os.remove(temp_outfile)
    return cobra_model


def find_outliers_iqr(df):
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    outlier_mask = (df < (Q1 - 1.5 * IQR)) | (df > (Q3 + 1.5 * IQR))
    return df.index[outlier_mask.any(axis=1)].tolist()


def set_publication_style(
    base_size: int = 10,
    base_family: str = 'sans-serif',
    width_mm: float = 90.0,
    height_ratio: float = 0.66,
):
    width_in = width_mm / 25.4
    height_in = width_in * height_ratio

    plt.style.use('default')
    sns.set_theme(style='white')

    plt.rcParams.update({
        # Fonts
        'font.size': base_size,
        'font.family': base_family,
        'axes.titlesize': base_size + 1,
        'axes.titleweight': 'medium',
        'axes.titlepad': 8,
        'axes.labelsize': base_size,
        'axes.labelweight': 'medium',
        'xtick.labelsize': max(base_size - 1, 6),
        'ytick.labelsize': max(base_size - 1, 6),
        'legend.fontsize': max(base_size - 1, 6),
        'legend.title_fontsize': base_size,

        # Figure sizing and backgrounds
        'figure.figsize': (width_in, height_in),
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',

        # Axes and spines (minimal)
        'axes.grid': False,
        'axes.edgecolor': 'black',
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.spines.left': True,
        'axes.spines.bottom': True,

        # Ticks (clear and outward)
        'xtick.bottom': True,
        'ytick.left': True,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.major.size': 3.5,
        'ytick.major.size': 3.5,
        'xtick.minor.visible': False,
        'ytick.minor.visible': False,
        'xtick.color': 'black',
        'ytick.color': 'black',

        # Legend (clean)
        'legend.frameon': False,

        # Titles
        'axes.titlelocation': 'center',

        # Layout and saving
        'figure.autolayout': False,
        'figure.constrained_layout.use': True,
        'savefig.bbox': 'tight',
        'savefig.dpi': 300,
        'savefig.facecolor': 'white',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })


def run_umap(cobra_model, output_dir, output_viz_dir, flux_comparison_file,
             mode='all', random_state=42):
    """Run UMAP dimensionality reduction.

    Parameters
    ----------
    cobra_model          : COBRApy Model (kept for API compatibility)
    output_dir           : directory containing flux1.csv / flux2.csv
    output_viz_dir       : output directory for UMAP results
    flux_comparison_file : CSV of reactions to include
    mode                 : 'all' or 'df'
    random_state         : RNG seed

    Returns
    -------
    target_reactions : list of reaction IDs used
    """
    label_str = 'All Flux' if mode == 'all' else 'DF reactions only'
    logging.info(f'\n[UMAP-{mode}] {label_str}')

    # Reaction selection: all common reactions from flux1 ∩ flux2
    df1_full = pd.read_csv(f'{output_dir}/flux1.csv', index_col=0)
    df2_full = pd.read_csv(f'{output_dir}/flux2.csv', index_col=0)
    common_all = list(df1_full.index.intersection(df2_full.index))
    logging.info(f'  Common reactions (flux1 ∩ flux2): {len(common_all)}')

    if mode == 'df':
        ref_df = pd.read_csv(flux_comparison_file, index_col=0)
        target_reactions = [r for r in common_all if r in ref_df.index]
        logging.info(f'  Target reactions after DF filter: {len(target_reactions)}')
    else:
        target_reactions = common_all
        logging.info(f'  Target reactions (all): {len(target_reactions)}')

    # Load & merge flux matrices
    df1 = df1_full.reindex(target_reactions).dropna(how='all')
    df2 = df2_full.reindex(target_reactions).dropna(how='all')
    common_rxns = df1.index.intersection(df2.index)
    df1 = df1.loc[common_rxns]
    df2 = df2.loc[common_rxns]

    c_cols = list(df1.columns)
    t_cols = list(df2.columns)
    merged = pd.concat([df1.T, df2.T], axis=0).fillna(0)

    rv     = merged.var(axis=0)
    merged = merged.loc[:, rv > 0]
    logging.info(f'  Samples: {len(c_cols)} Control + {len(t_cols)} Test  |  '
                 f'Reactions after zero-var filter: {merged.shape[1]}')

    conditions = np.array(['Control'] * len(c_cols) + ['Test'] * len(t_cols))

    os.makedirs(output_viz_dir, exist_ok=True)

    # StandardScaler + PCA + UMAP
    scaler    = StandardScaler()
    X_scaled  = scaler.fit_transform(merged.values)

    n_pcs     = min(30, merged.shape[0] - 1, merged.shape[1])
    pca_model = PCA(n_components=n_pcs, random_state=random_state)
    pca_scores = pca_model.fit_transform(X_scaled)

    reducer   = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=random_state)
    embedding = reducer.fit_transform(pca_scores)

    out_df = pd.DataFrame(embedding, index=merged.index, columns=['UMAP1', 'UMAP2'])
    out_df['condition'] = conditions
    out_df.to_csv(os.path.join(output_viz_dir, f'UMAP_{mode}_raw_data.csv'))

    # Save models for MOMA projection
    joblib.dump(scaler,    os.path.join(output_viz_dir, f'umap_{mode}_scaler.pkl'))
    joblib.dump(pca_model, os.path.join(output_viz_dir, f'umap_{mode}_pca.pkl'))
    joblib.dump(reducer,   os.path.join(output_viz_dir, f'umap_{mode}_model.pkl'))

    logging.info(f'  Saved: UMAP_{mode}_raw_data.csv | scaler/pca/umap pkl')
    return target_reactions


def run_r_umap(output_dir, output_viz_dir,
               df_file=None, mode='both', r_script=None,
               c1_color=None, c2_color=None,
               c1_label=None, c2_label=None):
    """Call flux_dimred.R via subprocess, then fit sklearn models from R output.

    R handles: PCA → UMAP → visualization → CSV/JSON outputs.
    Python reads R's merged flux matrix + optimal params, fits
    StandardScaler + PCA + UMAP with those exact parameters, and saves pkl
    files for MOMA projection.

    Parameters
    ----------
    output_dir        : directory containing flux1.csv / flux2.csv
    output_viz_dir    : where R writes results (and pkl files are saved)
    df_file           : Differential_fluxes.csv (optional, for df mode)
    mode              : 'all' | 'df' | 'both'
    r_script          : path to flux_dimred.R (auto-detected if None)

    Returns
    -------
    dict  {mode: target_reactions_list}
    """
    import subprocess

    if r_script is None:
        r_script = _r_script('flux_dimred.R')
    if not os.path.exists(r_script):
        raise FileNotFoundError(f'R script not found: {r_script}')

    os.makedirs(output_viz_dir, exist_ok=True)

    cmd = [
        _find_rscript(), r_script,
        '--output-dir',     output_dir,
        '--output-viz-dir', output_viz_dir,
        '--mode',           mode,
    ]
    if df_file and os.path.exists(df_file):
        cmd += ['--df-file', df_file]
    if c1_color:
        cmd += ['--c1-color', c1_color]
    if c2_color:
        cmd += ['--c2-color', c2_color]
    if c1_label:
        cmd += ['--c1-label', c1_label]
    if c2_label:
        cmd += ['--c2-label', c2_label]

    logging.info('[R] Running: %s', ' '.join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError('flux_dimred.R exited with code %d' % result.returncode)

    logging.info('[R] Done.')

    # ------------------------------------------------------------------
    # Fit sklearn models from R outputs (for MOMA projection)
    # ------------------------------------------------------------------
    merged_csv = os.path.join(output_viz_dir, 'merged_flux.csv')
    if not os.path.exists(merged_csv):
        raise FileNotFoundError(f'R did not produce: {merged_csv}')
    mat_corr = pd.read_csv(merged_csv, index_col=0)   # features × samples

    target_reactions = {}
    modes_to_fit = ['all', 'df'] if mode == 'both' else [mode]

    for m in modes_to_fit:
        params_json = os.path.join(output_viz_dir, f'umap_{m}_params.json')
        if not os.path.exists(params_json):
            logging.warning('Params file not found (R skipped mode=%s): %s', m, params_json)
            continue

        with open(params_json) as fp:
            params = json.load(fp)

        n_pcs    = int(params['n_pcs'])
        best_nn  = int(params['best_nn'])
        best_md  = float(params['best_dist'])
        features = params['features']

        target_reactions[m] = features

        # Subset + transpose: samples × features
        mat_m = mat_corr.reindex(features).fillna(0).T   # (n_samples, n_features)

        logging.info('[sklearn-%s] n_samples=%d  n_features=%d  n_pcs=%d  '
                     'n_neighbors=%d  min_dist=%s',
                     m, mat_m.shape[0], mat_m.shape[1], n_pcs, best_nn, best_md)

        scaler    = StandardScaler()
        X_scaled  = scaler.fit_transform(mat_m.values)

        pca_model = PCA(n_components=min(n_pcs, mat_m.shape[0] - 1, mat_m.shape[1]),
                        random_state=42)
        pca_scores = pca_model.fit_transform(X_scaled)

        reducer   = umap.UMAP(n_neighbors=best_nn, min_dist=best_md,
                              n_components=2, random_state=42)
        embedding = reducer.fit_transform(pca_scores)   # (n_samples, 2)

        joblib.dump(scaler,    os.path.join(output_viz_dir, f'umap_{m}_scaler.pkl'))
        joblib.dump(pca_model, os.path.join(output_viz_dir, f'umap_{m}_pca.pkl'))
        joblib.dump(reducer,   os.path.join(output_viz_dir, f'umap_{m}_model.pkl'))

        # Save Python-space coordinates as a SEPARATE reference file.
        # R's UMAP_*_raw_data.csv is kept intact so UMAP_*_.png (flux viz) is unchanged.
        # MOMA projection uses this Python-space reference for correct coordinate alignment.
        ref_csv = os.path.join(output_viz_dir, f'UMAP_{m}_raw_data.csv')
        python_ref_csv = os.path.join(output_viz_dir, f'UMAP_{m}_python_ref.csv')
        if os.path.exists(ref_csv):
            ref_df = pd.read_csv(ref_csv, index_col=0)
            for i, sample in enumerate(mat_m.index):
                if sample in ref_df.index:
                    ref_df.loc[sample, 'UMAP1'] = embedding[i, 0]
                    ref_df.loc[sample, 'UMAP2'] = embedding[i, 1]
            ref_df.to_csv(python_ref_csv)
            logging.info('[sklearn-%s] Saved Python-space ref: %s',
                         m, os.path.basename(python_ref_csv))

        logging.info('[sklearn-%s] Saved: scaler / pca / umap pkl', m)

    return target_reactions


def _confidence_ellipse(ax, x, y, color, level=0.80, zorder=0):
    """Normal-theory confidence ellipse, matching ggplot stat_ellipse(type='norm').

    Draws a faint filled ellipse plus a dashed outline (same style as the R
    UMAP/PCA panels in flux_dimred.R / redraw_umap_panels.R).
    """
    from matplotlib.patches import Ellipse
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 3:
        return
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    n_std = np.sqrt(-2.0 * np.log(1.0 - level))   # chi-square radius (df = 2)
    width, height = 2 * n_std * np.sqrt(vals)
    cx, cy = x.mean(), y.mean()
    ax.add_patch(Ellipse((cx, cy), width, height, angle=theta,
                         facecolor=color, edgecolor='none', alpha=0.06, zorder=zorder))
    ax.add_patch(Ellipse((cx, cy), width, height, angle=theta,
                         facecolor='none', edgecolor=color,
                         linestyle='--', linewidth=0.8, zorder=zorder + 1))


def visualize_flux_umap(output_viz_dir, mode='all', c1_color=None, c2_color=None,
                        c1_label=None, c2_label=None, ellipse_level=0.80):
    """Matplotlib UMAP scatter mirroring flux_dimred.R / redraw_umap_panels.R.

    Reads ``UMAP_{mode}_raw_data.csv`` and draws condition-colored points with
    80% normal-theory confidence ellipses (dashed outline + faint fill), a plain
    title, a square panel, and a right-side "Condition" legend. Writes
    ``UMAP_{mode}.png`` (+ .svg) — the same output as the R backend.
    """
    set_publication_style()
    os.makedirs(output_viz_dir, exist_ok=True)

    df = pd.read_csv(os.path.join(output_viz_dir, f'UMAP_{mode}_raw_data.csv'),
                     index_col=0)

    color_map = {
        'Control': c1_color or CONDITION_COLORS['Control'],
        'Test':    c2_color or CONDITION_COLORS['Test'],
    }

    # Legend labels: prefer override, then the R-written 'label' column, then condition
    def _disp(cond, override):
        if override:
            return override
        if 'label' in df.columns:
            vals = df.loc[df['condition'] == cond, 'label'].dropna().unique()
            if len(vals):
                return str(vals[0])
        return cond
    label_map = {'Control': _disp('Control', c1_label),
                 'Test':    _disp('Test', c2_label)}

    fig, ax = plt.subplots(figsize=(5, 5))
    for cond in ['Control', 'Test']:
        m = df['condition'] == cond
        if not m.any():
            continue
        x, y = df.loc[m, 'UMAP1'].values, df.loc[m, 'UMAP2'].values
        _confidence_ellipse(ax, x, y, color_map[cond], level=ellipse_level)
        ax.scatter(x, y, c=color_map[cond], marker='o', s=28, alpha=0.85,
                   linewidths=0, label=label_map[cond], zorder=3)

    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('UMAP')
    ax.set_box_aspect(1)
    ax.legend(title='Condition', frameon=False,
              loc='center left', bbox_to_anchor=(1.02, 0.5))

    fig.savefig(os.path.join(output_viz_dir, f'UMAP_{mode}.svg'),
                format='svg', bbox_inches='tight')
    fig.savefig(os.path.join(output_viz_dir, f'UMAP_{mode}.png'),
                format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info('Saved: UMAP_%s.svg / .png', mode)


def run_moma_results_umap(output_dir, output_viz_dir, target_reactions,
                          mode='df', use_r=True, moma_subdir=''):
    """Project MOMA targeting results into UMAP space.

    Parameters
    ----------
    output_dir       : directory containing flux2.csv and targeting_results/
    output_viz_dir   : directory containing UMAP models and ComBat CSV
    target_reactions : list of reaction IDs (used by Python path;
                       R reads features from umap_{mode}_params.json)
    mode             : UMAP reference mode ('df' or 'all')
    use_r            : if True, call moma_project.R — coordinates match
                         UMAP_{mode}_raw_data.csv (R/ggplot2 space)
                       if False, use Python sklearn models — coordinates
                         match UMAP_{mode}_python_ref.csv (sklearn space)
    moma_subdir      : subdir under targeting_results/ to read MOMA files,
                       and subdir under output_viz_dir/ to write results.
                       Default '' = root level (original behavior).
    """
    if use_r:
        import subprocess
        r_script = _r_script('moma_project.R')
        if not os.path.exists(r_script):
            raise FileNotFoundError(f'moma_project.R not found: {r_script}')
        cmd = [
            _find_rscript(), r_script,
            '--output-dir',     output_dir,
            '--output-viz-dir', output_viz_dir,
            '--mode',           mode,
        ]
        if moma_subdir:
            cmd += ['--moma-subdir', moma_subdir]
        logging.info('[MOMA-R] Running: %s', ' '.join(cmd))
        result = subprocess.run(cmd)
        if result.returncode != 0:
            raise RuntimeError('moma_project.R exited with code %d' % result.returncode)
        logging.info('[MOMA-R] Done.')
        return

    # ---- Python (sklearn) path -----------------------------------------------
    import joblib

    # Load sklearn models
    scaler    = joblib.load(os.path.join(output_viz_dir, f'umap_{mode}_scaler.pkl'))
    pca_model = joblib.load(os.path.join(output_viz_dir, f'umap_{mode}_pca.pkl'))
    reducer   = joblib.load(os.path.join(output_viz_dir, f'umap_{mode}_model.pkl'))

    # Load features from params JSON
    params_json = os.path.join(output_viz_dir, f'umap_{mode}_params.json')
    with open(params_json) as fp:
        params = json.load(fp)
    features = params['features']

    # Control center from Python-space reference
    python_ref_csv = os.path.join(output_viz_dir, f'UMAP_{mode}_python_ref.csv')
    ref_df   = pd.read_csv(python_ref_csv, index_col=0)
    ctrl_df  = ref_df[ref_df['condition'] == 'Control']
    center_x = np.median(ctrl_df['UMAP1'].values)
    center_y = np.median(ctrl_df['UMAP2'].values)

    # Load merged flux baseline and raw flux2 for delta correction
    merged_df = pd.read_csv(os.path.join(output_viz_dir, 'merged_flux.csv'),
                            index_col=0)
    flux2_df  = pd.read_csv(os.path.join(output_dir, 'flux2.csv'),
                            index_col=0).fillna(0)

    # Project each MOMA result file
    moma_dir   = os.path.join(output_dir, 'targeting_results')
    moma_files = glob.glob(os.path.join(moma_dir, 'MOMA_target_results_*.csv'))
    logging.info('[MOMA-Python-%s] Projecting %d files ...', mode, len(moma_files))

    for moma_file in moma_files:
        base_noext  = os.path.splitext(os.path.basename(moma_file))[0]
        sample_name = base_noext.replace('MOMA_target_results_', '')

        df_moma  = pd.read_csv(moma_file, index_col=0)   # reactions × genes
        df_feat  = df_moma.reindex(features).fillna(0)   # features × genes
        moma_mat = df_feat.T.values                       # genes × features

        # Delta correction
        if sample_name in merged_df.columns and sample_name in flux2_df.columns:
            wt_combat = merged_df.reindex(features)[sample_name].fillna(0).values
            wt_raw    = flux2_df.reindex(features)[sample_name].fillna(0).values
            X_input   = (moma_mat - wt_raw[np.newaxis, :]) + wt_combat[np.newaxis, :]
        else:
            logging.warning('[MOMA-Python] %s not in combat/flux2 — raw projection',
                            sample_name)
            X_input = moma_mat

        # Scaler → PCA → UMAP
        X_scaled    = scaler.transform(X_input)
        pca_scores  = pca_model.transform(X_scaled)
        umap_coords = reducer.transform(pca_scores)   # genes × 2

        gene_names = list(df_moma.columns)
        result_df  = pd.DataFrame({
            'X': umap_coords[:, 0],
            'Y': umap_coords[:, 1],
            'distance_to_center': np.sqrt((umap_coords[:, 0] - center_x) ** 2 +
                                          (umap_coords[:, 1] - center_y) ** 2),
        }, index=gene_names)

        out_csv = os.path.join(output_viz_dir,
                               f'UMAP_Results_{base_noext}.csv')
        result_df.to_csv(out_csv)

    logging.info('[MOMA-Python-%s] Done.', mode)


def visualize_moma_umap(output_viz_dir, mode='df', use_r=False):
    """Visualize MOMA simulation results.

    Parameters
    ----------
    output_viz_dir : directory containing UMAP_Results_*.csv
    mode           : UMAP reference mode ('df' or 'all')
    use_r          : if True, delegate to moma_viz.R (ggplot2 + ggrepel)
    """
    if use_r:
        import subprocess
        r_script = _r_script('moma_viz.R')
        cmd = [_find_rscript(), r_script,
               '--output-viz-dir', output_viz_dir,
               '--mode',           mode]
        # gene-id → HGNC symbol map (run-local Human ENSG map if present, else bundled Entrez)
        gene_map = _effective_gene_map(output_viz_dir)
        if os.path.exists(gene_map):
            cmd += ['--gene-map', gene_map]
        logging.info('[R] Running: %s', ' '.join(cmd))
        result = subprocess.run(cmd)
        if result.returncode != 0:
            raise RuntimeError('moma_viz.R exited with code %d' % result.returncode)
        logging.info('[R] MOMA visualization done.')
        return

    # ---- Python (matplotlib) path — mirrors moma_viz.R style --------------
    import io, sys
    import matplotlib.patheffects as pe
    try:
        from adjustText import adjust_text          # optional label de-overlap
    except ImportError:
        adjust_text = None

    set_publication_style()

    CTRL_COLOR   = '#009E73'   # Control (Wong bluish green)
    TARGET_COLOR = '#E69F00'   # Top-10 targets (orange)
    SAMPLE_COLOR = '#D55E00'   # Sample / toward-Control arrows (vermilion)
    CENTER_COLOR = '#333333'   # Control Center cross
    AWAY_COLOR   = '#B0B0B0'   # arrows away from Control / background MOMA points

    # gene-id → HGNC symbol map (run-local Human ENSG map if present, else bundled Entrez)
    gene_map_file = _effective_gene_map(output_viz_dir)
    entrez2symbol = {}
    if os.path.exists(gene_map_file):
        gm = pd.read_csv(gene_map_file, dtype=str)
        if {'entrez', 'symbol'}.issubset(gm.columns):
            entrez2symbol = dict(zip(gm['entrez'], gm['symbol']))

    def _to_symbol(g):
        s = entrez2symbol.get(str(g))
        return s if isinstance(s, str) and s else str(g)

    # Reference UMAP (python_ref preferred, else R raw data)
    python_ref_csv = os.path.join(output_viz_dir, f'UMAP_{mode}_python_ref.csv')
    raw_ref_csv    = os.path.join(output_viz_dir, f'UMAP_{mode}_raw_data.csv')
    ref_csv        = python_ref_csv if os.path.exists(python_ref_csv) else raw_ref_csv
    umap_ref_df    = pd.read_csv(ref_csv, index_col=0)

    ctrl_mask = umap_ref_df['condition'].values == 'Control'
    ctrl_x    = umap_ref_df.loc[ctrl_mask, 'UMAP1'].values
    ctrl_y    = umap_ref_df.loc[ctrl_mask, 'UMAP2'].values
    center_x  = np.median(ctrl_x)
    center_y  = np.median(ctrl_y)

    result_files = glob.glob(os.path.join(output_viz_dir, 'UMAP_Results_*.csv'))
    target_rows = []

    for result_file in result_files:
        basename    = os.path.basename(result_file).replace('UMAP_Results_', '').split('.')[0]
        sample_name = basename.replace('MOMA_target_results_', '')
        tmp_df      = pd.read_csv(result_file, index_col=0)

        # Rank top-10 by recovery toward Control (recovery_proj_pca desc if present,
        # else fall back to distance_to_center asc)
        if 'recovery_proj_pca' in tmp_df.columns:
            top = tmp_df.sort_values('recovery_proj_pca', ascending=False).head(10)
        else:
            top = tmp_df.sort_values('distance_to_center', ascending=True).head(10)

        sample_xy = None
        if sample_name in umap_ref_df.index:
            sample_xy = (umap_ref_df.loc[sample_name, 'UMAP1'],
                         umap_ref_df.loc[sample_name, 'UMAP2'])

        fig, ax = plt.subplots(figsize=(7.2, 6.2))

        _confidence_ellipse(ax, ctrl_x, ctrl_y, CTRL_COLOR, level=0.80)
        ax.scatter(ctrl_x, ctrl_y, color=CTRL_COLOR, s=18, alpha=0.7,
                   linewidths=0, label='Control', zorder=2)
        # All MOMA gene KO results (gray cloud)
        ax.scatter(tmp_df['X'], tmp_df['Y'], color=AWAY_COLOR, s=10, alpha=0.3,
                   linewidths=0, zorder=1)

        # Arrows: sample → KO position for top-10 (vermilion = toward Control)
        if sample_xy is not None:
            sx, sy = sample_xy
            d_sample = np.hypot(sx - center_x, sy - center_y)
            for gid, r in top.iterrows():
                if 'recovery_proj_pca' in top.columns:
                    toward = r['recovery_proj_pca'] > 0
                else:
                    toward = np.hypot(r['X'] - center_x, r['Y'] - center_y) < d_sample
                ax.annotate('', xy=(r['X'], r['Y']), xytext=(sx, sy),
                            arrowprops=dict(arrowstyle='->',
                                            color=SAMPLE_COLOR if toward else AWAY_COLOR,
                                            alpha=0.55 if toward else 0.4, lw=0.6),
                            zorder=3)

        ax.scatter(top['X'], top['Y'], color=TARGET_COLOR, s=40, alpha=0.9,
                   linewidths=0, label='Top-10 targets', zorder=4)
        # Cross markers: Control Center + Sample
        ax.scatter(center_x, center_y, color=CENTER_COLOR, marker='+',
                   s=120, linewidths=1.6, zorder=6, label='Control Center')
        if sample_xy is not None:
            ax.scatter(sample_xy[0], sample_xy[1], color=SAMPLE_COLOR, marker='+',
                       s=120, linewidths=1.6, zorder=6, label=f'Sample({sample_name})')

        # Gene labels (HGNC symbols, italic, white halo + leader lines via adjustText)
        texts = [ax.text(r['X'], r['Y'], _to_symbol(gid), fontsize=9, fontstyle='italic',
                         path_effects=[pe.withStroke(linewidth=2, foreground='white')])
                 for gid, r in top.iterrows()]

        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_title('MOMA Simulation results')
        ax.set_box_aspect(1)
        ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1.02, 0.5))

        if adjust_text is not None:
            _old_stdout = sys.stdout
            sys.stdout = io.StringIO()
            try:
                adjust_text(texts, ax=ax,
                            arrowprops=dict(arrowstyle='-', color='#8C8C8C', lw=0.4))
            finally:
                sys.stdout = _old_stdout

        fig.savefig(os.path.join(output_viz_dir, f'UMAP_{basename}.svg'),
                    format='svg', bbox_inches='tight')
        fig.savefig(os.path.join(output_viz_dir, f'UMAP_{basename}.png'),
                    format='png', dpi=300, bbox_inches='tight')
        plt.close(fig)

        for gid in top.index:
            target_rows.append((str(gid), _to_symbol(gid), sample_name))

    pd.DataFrame(target_rows, columns=['Gene', 'Symbol', 'Sample']).to_csv(
        os.path.join(output_viz_dir, 'target_genes.csv'), index=False)
    return


def _l2_score_from_moma(output_dir, sample_name, gene_index):
    """Fallback per-knockout score: L2 displacement from the patient's WT flux.

    Used only when ``recovery_proj_pca`` is unavailable. For each single-gene
    knockout, computes the Euclidean (L2) distance of its perturbed flux vector
    (``targeting_results/MOMA_target_results_<id>.csv``, reactions × genes; MOMA
    already restricts to non-boundary reactions) from the patient's wild-type
    flux (``flux2/Flux_prediction_<id>.csv``). Larger displacement = stronger
    metabolic perturbation. Returns a Series of scores aligned to ``gene_index``,
    or None if the required files are missing.
    """
    moma_file = os.path.join(output_dir, 'targeting_results',
                             'MOMA_target_results_%s.csv' % sample_name)
    wt_files = glob.glob(os.path.join(output_dir, 'flux2',
                                      'Flux_prediction_%s.csv' % sample_name))
    if not os.path.exists(moma_file) or not wt_files:
        return None
    moma = pd.read_csv(moma_file, index_col=0)                  # reactions × genes
    wt = pd.read_csv(wt_files[0], index_col=0, header=None).iloc[:, 0]
    wt_vec = wt.reindex(moma.index).fillna(0.0).values.reshape(-1, 1)
    ko = moma.fillna(0.0).values                                # reactions × genes
    l2 = np.sqrt(((ko - wt_vec) ** 2).sum(axis=0))              # per-gene L2
    s = pd.Series(l2, index=[str(c) for c in moma.columns])
    return s.reindex([str(g) for g in gene_index])


def export_moma_gene_scores_full(output_dir, output_viz_dir, mode='df'):
    """Persist the FULL per-gene-per-patient MOMA target score table.

    Writes ``viz/moma_gene_scores_full.csv`` (columns: Gene, Symbol, Sample,
    score, rank) — the complete per-patient ranking required for cross-cohort
    rank aggregation (``cohort_target_ranking``). The existing top-10
    ``target_genes.csv`` is left unchanged.

    ``score`` definition
    --------------------
    Per patient, every single-gene knockout is scored by how much it moves the
    perturbed flux state toward the Control center (higher = better target):

      * Preferred — ``recovery_proj_pca`` from ``viz/UMAP_Results_*.csv`` (linear
        projection onto the patient→Control axis in PCA space; > 0 = toward
        Control).
      * Fallback — when ``recovery_proj_pca`` is absent, the L2 displacement of
        each knockout's flux vector from the patient's wild-type flux over
        non-boundary reactions (see ``_l2_score_from_moma``).

    ``rank`` is 1 = best (highest score) within each patient. Reads the
    per-patient ``UMAP_Results_*.csv`` files, so it is independent of whether the
    R or matplotlib visualization backend was used. Returns the output CSV path,
    or None if no per-patient MOMA result files are found.
    """
    gene_map_file = _effective_gene_map(output_viz_dir)
    entrez2symbol = {}
    if os.path.exists(gene_map_file):
        gm = pd.read_csv(gene_map_file, dtype=str)
        if {'entrez', 'symbol'}.issubset(gm.columns):
            entrez2symbol = dict(zip(gm['entrez'], gm['symbol']))

    def _to_symbol(g):
        s = entrez2symbol.get(str(g))
        return s if isinstance(s, str) and s else str(g)

    result_files = sorted(glob.glob(os.path.join(output_viz_dir, 'UMAP_Results_*.csv')))
    parts = []
    for result_file in result_files:
        base_noext  = os.path.basename(result_file).replace('UMAP_Results_', '').split('.csv')[0]
        sample_name = base_noext.replace('MOMA_target_results_', '')
        df_res = pd.read_csv(result_file, index_col=0)

        if 'recovery_proj_pca' in df_res.columns and df_res['recovery_proj_pca'].notna().any():
            score = df_res['recovery_proj_pca'].astype(float).values
        else:
            s = _l2_score_from_moma(output_dir, sample_name, df_res.index)
            if s is None:
                logging.warning('No recovery/MOMA score available for %s — skipping.', sample_name)
                continue
            score = s.values

        sdf = pd.DataFrame({
            'Gene':   [str(g) for g in df_res.index],
            'Symbol': [_to_symbol(g) for g in df_res.index],
            'Sample': sample_name,
            'score':  score,
        }).dropna(subset=['score'])
        if sdf.empty:
            continue
        sdf['rank'] = sdf['score'].rank(ascending=False, method='min').astype(int)
        parts.append(sdf.sort_values('rank'))

    if not parts:
        logging.warning('No per-patient MOMA score files found — '
                        'skipping moma_gene_scores_full.csv')
        return None

    full = pd.concat(parts, ignore_index=True)
    out_csv = os.path.join(output_viz_dir, 'moma_gene_scores_full.csv')
    full.to_csv(out_csv, index=False)
    logging.info('Saved: %s (%d gene-patient rows, %d patients)',
                 out_csv, len(full), full['Sample'].nunique())
    return out_csv


def annotate_de_flux_membership(output_dir, cobra_model, de_flux_reaction_ids):
    """Add/refresh ``de_flux`` and ``subsystem`` columns on cohort_target_ranking.csv,
    both from the model GPR (authoritative, so CONSISTENT within a GPR group — coupled
    isozymes share reactions): ``de_flux`` = any of the gene's reactions is in the
    differential-flux set; ``subsystem`` = the gene's most frequent reaction subsystem
    (its metabolic pathway). Idempotent; no-op if the ranking file is absent."""
    from collections import Counter
    ctr_csv = os.path.join(output_dir, 'cohort_target_ranking.csv')
    if not os.path.exists(ctr_csv) or os.path.getsize(ctr_csv) == 0:
        return None
    ctr = pd.read_csv(ctr_csv)
    de = set(str(r) for r in de_flux_reaction_ids)
    model_gene_ids = set(g.id for g in cobra_model.genes)
    # Enriched pathways with their best (lowest) enrichment FDR. When a gene spans several
    # subsystems, label it by its MOST-enriched pathway — the disease-relevant one — rather
    # than the subsystem with the most (often nucleotide-promiscuous) reactions. This makes
    # e.g. pyruvate kinase read "Glycolysis" (the dominant signal) not "Pyrimidine".
    enr_fdr = {}
    for fn in ('Up_regulated_pathways.csv', 'Down_regulated_pathways.csv'):
        p = os.path.join(output_dir, fn)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            edf = pd.read_csv(p, index_col=0)
            fcol = next((c for c in edf.columns if 'adjust' in c.lower() or c.lower() == 'fdr'), None)
            for path in edf.index.astype(str):
                fdr = float(edf.loc[path, fcol]) if fcol else 0.0
                enr_fdr[path] = min(enr_fdr.get(path, 1.0), fdr)

    def _rxns(gid):
        gid = str(gid)
        return list(cobra_model.genes.get_by_id(gid).reactions) if gid in model_gene_ids else []

    def _is_de(gid):
        return any(r.id in de for r in _rxns(gid))

    # generic catch-all subsystems: informative only when a gene has nothing more specific
    GENERIC = {'transport reactions', 'exchange/demand reactions', 'exchange reactions',
               'isolated', 'miscellaneous', 'pool reactions', 'artificial reactions', ''}

    def _subsystem(gid):
        subs = [r.subsystem.split(';')[0].strip() for r in _rxns(gid)
                if getattr(r, 'subsystem', '')]
        if not subs:
            return ''
        cnt = Counter(subs)
        enr = [s for s in cnt if s in enr_fdr]                       # enriched subsystems
        enr_spec = [s for s in enr if s.lower() not in GENERIC]      # ...that are specific
        if enr_spec:
            return min(enr_spec, key=lambda s: enr_fdr[s])          # most-enriched specific
        if enr:
            return min(enr, key=lambda s: enr_fdr[s])               # enriched (generic ok)
        spec = [s for s in cnt if s.lower() not in GENERIC]         # no enrichment: prefer
        pool = spec if spec else list(cnt)                          # specific by frequency
        return max(pool, key=lambda s: cnt[s])

    def _subsystem_frequent(gid):
        subs = [r.subsystem.split(';')[0].strip() for r in _rxns(gid)
                if getattr(r, 'subsystem', '')]
        return Counter(subs).most_common(1)[0][0] if subs else ''

    genes = ctr['Gene'].astype(str)
    ctr['de_flux']            = [bool(_is_de(g)) for g in genes]
    ctr['subsystem']          = [_subsystem(g) for g in genes]          # most-enriched pathway
    ctr['subsystem_frequent'] = [_subsystem_frequent(g) for g in genes]  # most-frequent subsystem
    ctr.to_csv(ctr_csv, index=False)
    logging.info('Annotated de_flux + subsystem(enriched/frequent) on %s (%d/%d in DE flux).',
                 ctr_csv, int(ctr['de_flux'].sum()), len(ctr))
    return ctr_csv


def visualize_target_oncoprint(output_dir, output_viz_dir, top_k=10,
                               max_genes=40, fallback_n=20):
    """Oncoprint-style matrix of the cohort target ranking → viz/target_oncoprint.*

    Rows = target genes, columns = patients. A cell is filled when the gene is
    among that patient's top-k knockouts (within-patient ``rank <= top_k``), and
    shaded by that rank with a colorblind-safe sequential map (viridis; darkest =
    rank 1). Non-top-k / absent cells are left light grey. Marginal bars show
    per-gene recurrence (right) and per-patient hit count (top); GPR-coupled genes
    (shared ``gpr_group``) are placed adjacently and bracketed so coupled subunits
    read as one unit. Writes ``target_oncoprint.png`` (300 dpi) and ``.svg`` into
    ``output_viz_dir`` — matching the other figures. Deterministic (stable sorts).

    Inputs (not modified): ``<output_dir>/cohort_target_ranking.csv`` and
    ``<output_viz_dir>/moma_gene_scores_full.csv``. Gene selection: all genes with
    ``significant == True``; if none, fall back to the top ``fallback_n`` by
    ``rra_score`` with an explanatory subtitle; if more than ``max_genes`` are
    significant, cap at the top ``max_genes`` by ``rra_score``. Logs a warning and
    returns (no crash) when an input file is missing or empty.

    Parameters
    ----------
    output_dir     : directory holding cohort_target_ranking.csv
    output_viz_dir : viz directory holding moma_gene_scores_full.csv (and output)
    top_k          : within-patient top-k cutoff for a filled cell (default 10)
    max_genes      : cap on significant genes shown (default 40)
    fallback_n     : number of genes shown when none are significant (default 20)
    """
    from matplotlib.colors import Normalize
    from matplotlib import cm
    from matplotlib.lines import Line2D

    ctr_csv  = os.path.join(output_dir, 'cohort_target_ranking.csv')
    full_csv = os.path.join(output_viz_dir, 'moma_gene_scores_full.csv')
    for f in (ctr_csv, full_csv):
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            logging.warning('Oncoprint skipped — missing/empty input: %s', f)
            return None

    ctr  = pd.read_csv(ctr_csv)
    full = pd.read_csv(full_csv)
    if ctr.empty or full.empty:
        logging.warning('Oncoprint skipped — empty cohort ranking or score table.')
        return None

    ctr['Gene']  = ctr['Gene'].astype(str)
    full['Gene'] = full['Gene'].astype(str)
    k = max(int(top_k), 1)

    # ---- Gene selection (significant; else top-N fallback; cap at max_genes) ---
    sig_mask = ctr['significant'].astype(str).str.lower().isin(['true', '1', '1.0'])
    fallback = False
    if sig_mask.any():
        sel = ctr[sig_mask].copy()
        if len(sel) > max_genes:
            sel = sel.sort_values('rra_score', kind='mergesort').head(max_genes)
    else:
        fallback = True
        sel = ctr.sort_values('rra_score', kind='mergesort').head(fallback_n).copy()

    # ---- Row order: recurrence_prop desc, rra_score asc; GPR groups adjacent ---
    sel['gpr_group'] = sel['gpr_group'].fillna('').astype(str)
    sel['_grp'] = [g.strip() if g.strip() else '__solo_%s' % gene
                   for g, gene in zip(sel['gpr_group'], sel['Gene'])]
    sel = sel.sort_values(['recurrence_prop', 'rra_score'],
                          ascending=[False, True], kind='mergesort').reset_index(drop=True)
    sel['_pos'] = np.arange(len(sel))
    sel['_grp_first'] = sel.groupby('_grp')['_pos'].transform('min')
    sel = sel.sort_values(['_grp_first', '_pos'], kind='mergesort').reset_index(drop=True)

    row_genes = list(sel['Gene'])
    # Row labels: prefer cohort Symbol, else gene id
    def _label(row):
        s = row.get('Symbol')
        return str(s) if isinstance(s, str) and s and s.lower() != 'nan' else str(row['Gene'])
    row_labels  = [_label(r) for _, r in sel.iterrows()]
    row_recur   = sel['recurrence_prop'].astype(float).values
    row_sig     = sel['significant'].astype(str).str.lower().isin(['true', '1', '1.0']).values
    # DE-flux membership (gene's reactions are differentially active): drives the '*'.
    row_deflux  = (sel['de_flux'].astype(str).str.lower().isin(['true', '1', '1.0']).values
                   if 'de_flux' in sel.columns else np.zeros(len(sel), dtype=bool))

    # ---- Rank matrix (genes x patients) ---------------------------------------
    fsel = full[full['Gene'].isin(set(row_genes))]
    rankmat = fsel.pivot_table(index='Gene', columns='Sample', values='rank', aggfunc='min')
    patients_all = list(full['Sample'].astype(str).unique())
    rankmat = rankmat.reindex(index=row_genes, columns=patients_all)

    M = rankmat.values.astype(float)
    M[~np.isfinite(M)] = np.nan
    M[M > k] = np.nan                         # only top-k cells are filled

    # ---- Column order: patient hit count desc (stable) ------------------------
    col_hits_all = np.isfinite(M).sum(axis=0)
    order = np.argsort(-col_hits_all, kind='stable')
    patients = [patients_all[i] for i in order]
    M = M[:, order]
    col_hits = np.isfinite(M).sum(axis=0)
    n, m = M.shape

    import copy as _copy
    masked = np.ma.masked_invalid(M)
    cmap = _copy.copy(cm.get_cmap('viridis'))
    cmap.set_bad('#f0f0f0')                   # non-top-k / absent = light grey
    norm = Normalize(vmin=1, vmax=max(k, 2))

    # ---- Layout (fixed fractional axes; figure size scales with matrix) -------
    set_publication_style()
    fig_w = max(8.0, 0.24 * m + 4.8)
    fig_h = max(4.5, 0.32 * n + 2.6)
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')
    fig.set_constrained_layout(False)         # manual add_axes layout
    ax_left  = fig.add_axes([0.015, 0.085, 0.150, 0.70])   # labels + GPR brackets
    ax_main  = fig.add_axes([0.175, 0.085, 0.620, 0.70])   # heatmap
    ax_top   = fig.add_axes([0.175, 0.800, 0.620, 0.130], sharex=ax_main)  # patient bars
    ax_right = fig.add_axes([0.800, 0.085, 0.140, 0.70], sharey=ax_main)   # recurrence bars
    ax_cbar  = fig.add_axes([0.800, 0.840, 0.140, 0.028])  # rank colorbar

    # Heatmap
    im = ax_main.imshow(masked, aspect='auto', cmap=cmap, norm=norm,
                        origin='upper', interpolation='nearest')
    ax_main.set_xticks(np.arange(m))
    ax_main.set_xticklabels(patients, rotation=90, fontsize=max(6, 9 - m // 12))
    ax_main.set_yticks([])
    ax_main.set_xlim(-0.5, m - 0.5)
    ax_main.set_ylim(n - 0.5, -0.5)
    # thin cell separators
    ax_main.set_xticks(np.arange(-0.5, m, 1), minor=True)
    ax_main.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax_main.grid(which='minor', color='white', linewidth=0.6)
    ax_main.tick_params(which='minor', length=0)

    # Left: gene labels + GPR brackets (shares y with heatmap)
    ax_left.set_xlim(0, 1)
    ax_left.set_ylim(n - 0.5, -0.5)
    ax_left.axis('off')
    # '*' now marks DE-flux genes (their reactions are differentially active); all rows
    # shown are already FDR + effect-size significant, so a significance '*' was redundant.
    for i, (lab, de) in enumerate(zip(row_labels, row_deflux)):
        txt = lab + (' *' if de else '')
        ax_left.text(0.96, i, txt, ha='right', va='center',
                     fontsize=max(6, 9 - n // 18),
                     fontweight='bold' if de else 'normal',
                     color='black' if de else '#555555', fontstyle='italic')
    # GPR brackets for multi-member groups (contiguous rows by construction). The bracket
    # is labelled with the group's Human-GEM PATHWAY only (e.g. "Glycolysis"); the member
    # gene symbols are already the row labels, so repeating them on the bracket is redundant
    # and unreadable. Falls back to the family name if the pathway is unknown.
    def _group_name(symbols):
        pref = os.path.commonprefix(symbols)
        pref = ''.join(ch for ch in pref if ch.isalpha())
        return pref if len(pref) >= 3 else '/'.join(symbols)

    def _short_pathway(s):
        s = str(s).split(' / ')[0].strip()               # primary of "A / B" subsystems
        return '' if s.lower() in ('', 'nan') else (s if len(s) <= 22 else s[:21] + '…')
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
            xb = 0.17                                       # bracket just right of the label
            ax_left.plot([xb, xb], [i - 0.42, j + 0.42], color='#444444',
                         lw=1.1, clip_on=False)
            ax_left.plot([xb, xb + 0.045], [i - 0.42, i - 0.42], color='#444444',
                         lw=1.1, clip_on=False)
            ax_left.plot([xb, xb + 0.045], [j + 0.42, j + 0.42], color='#444444',
                         lw=1.1, clip_on=False)
            enr = _short_pathway(grp_sub_enr[i])          # most-enriched pathway
            frq = _short_pathway(grp_sub_frq[i])          # most-frequent subsystem
            lines = []
            if enr:
                lines.append('– %s –' % enr)     # en-dashes: most enriched
            if frq:
                lines.append('(%s)' % frq)                 # parentheses: most frequent
            glabel = '\n'.join(lines) if lines else _group_name(row_labels[i:j + 1])
            # label anchored to the LEFT, bracket to its right → a clear gap between them
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

    # Rank colorbar
    cb = fig.colorbar(im, cax=ax_cbar, orientation='horizontal')
    cb.set_label('within-patient rank (1 = top)', fontsize=7)
    cb.ax.tick_params(labelsize=6)
    cb.set_ticks([1, max(k, 2)])

    # Titles / legend
    fig.suptitle('Cohort target oncoprint', x=0.5, y=0.985, fontsize=13,
                 fontweight='medium')
    n_de = int(np.count_nonzero(row_deflux))
    sub = 'no significant targets — showing top %d by RRA' % fallback_n if fallback \
        else ('%d significant target genes (all FDR + effect-size significant)  '
              '•  * = gene is DE-flux included (%d/%d)' % (n, n_de, n))
    fig.text(0.175, 0.945, sub, fontsize=8.5, color='#333333')
    legend_handles = [Line2D([0], [0], marker='|', color='#444444', lw=0,
                             markersize=10, markeredgewidth=1.2, label='GPR-coupled group'),
                      Line2D([0], [0], marker='$*$', color='black', lw=0,
                             markersize=9, label='* = gene in DE flux (reactions differentially active)'),
                      Line2D([0], [0], color='none',
                             label='– pathway –  =  most enriched pathway'),
                      Line2D([0], [0], color='none',
                             label='(pathway)  =  most frequent subsystem')]
    # place the legend BELOW the rotated patient (column) labels so it never overlaps them
    # (labels hang under the heatmap; anchor the legend top just below the figure bottom)
    fig.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(0.015, -0.02),
               frameon=False, fontsize=7)

    out_png = os.path.join(output_viz_dir, 'target_oncoprint.png')
    out_svg = os.path.join(output_viz_dir, 'target_oncoprint.svg')
    fig.savefig(out_png, format='png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_svg, format='svg', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    logging.info('Saved: %s (+ .svg)  [%d genes x %d patients, %s]',
                 out_png, n, m, 'fallback top-%d' % fallback_n if fallback else 'significant')
    return out_png


def visualize_pathway_flux_heatmap(output_dir, output_viz_dir, level='reaction',
                                   pathways=None, de_only=False, top_n_per_pathway=5,
                                   fdr_cutoff=0.1):
    """Clustered flux heatmap of enriched-pathway reactions/genes → viz/pathway_flux_heatmap_*.

    ``top_n_per_pathway`` (default 5; None = no cap) keeps only the N most
    differentially-active features per pathway — ranked by ``|Effect size|`` (the
    DE-flux rank-biserial/Cohen's d, ties broken by ascending P) — so the figure
    stays compact and labeled instead of growing to hundreds/thousands of rows
    (e.g. a cohort with many enriched pathways). For the gene level a gene is
    ranked within a pathway by its strongest reaction effect.

    Reads existing pipeline outputs only (NO new analysis): the per-sample fluxes
    and DE-flux annotations in ``gene_subsystem_flux.csv`` and the enriched
    subsystems in ``Up_regulated_pathways.csv`` / ``Down_regulated_pathways.csv``.

    ``de_only`` selects WHICH reactions are shown within the enriched pathways:
      * ``False`` (default) — ALL flux-carrying reactions in the enriched
        subsystems. The enrichment selects the PATHWAYS; every reaction in them is
        plotted (so the ``Direction``/Status annotation reflects each reaction's
        mean-change direction, not its DE significance). Output filenames have no
        suffix.
      * ``True`` — restrict further to the differentially-active reactions only,
        i.e. those present in ``Differential_fluxes.csv`` (the DE-flux set). For
        the gene level the per-gene mean is then taken over the gene's DE reactions
        only. Output filenames get a ``_de`` suffix
        (``pathway_flux_heatmap_<level>_de.png``).

    Shared core for both levels:
      * Rows are restricted to features whose ``Subsystem`` is in the enriched set
        (or in ``pathways`` if given), then optionally to the DE-flux set
        (``de_only``).
      * Each row is z-scored across samples (rows with zero variance are dropped);
        the z-score color scale is a diverging, colorblind-safe map clipped to
        [-4, 4] and centered at 0.
      * Rows are hierarchically clustered (seaborn clustermap, deterministic
        Euclidean/average linkage) and annotated by ``Pathway`` (Subsystem) and
        ``Direction`` (Status: UP/DOWN/Reversed). Columns are ordered and SPLIT by
        ``Condition`` — Control (``flux1.csv`` samples) then Test (``flux2.csv``) —
        with a Condition color bar, so condition groups are visually separated
        (columns are not clustered, to keep the split intact).

    level='reaction'
        Rows = individual reactions (``Reaction`` column, de-duplicated); values are
        the per-sample reaction fluxes AS STORED (no aggregation). Row label =
        Reaction id. → ``viz/pathway_flux_heatmap_reaction.png`` (+ .svg).

    level='gene'
        Rows = genes; the reaction-level flux is aggregated to one value per gene
        per sample by the MEAN over that gene's reactions (within the enriched set).
        A gene appears once (GPR-coupled genes are kept as separate rows). Row
        annotation uses the gene's most frequent Subsystem/Status. Row label = HGNC
        Symbol (fallback to the gene id). → ``viz/pathway_flux_heatmap_gene.png``.

    Returns the PNG path, or None (with a logged warning) if a required input is
    missing or no usable rows remain. Deterministic.
    """
    import seaborn as sns
    import textwrap
    from matplotlib.patches import Patch

    level = level.lower()
    if level not in ('reaction', 'gene'):
        raise ValueError("level must be 'reaction' or 'gene'")

    gsf_csv   = os.path.join(output_dir, 'gene_subsystem_flux.csv')
    flux1_csv = os.path.join(output_dir, 'flux1.csv')
    flux2_csv = os.path.join(output_dir, 'flux2.csv')
    for f in (gsf_csv, flux1_csv, flux2_csv):
        if not os.path.exists(f):
            logging.warning('Pathway flux heatmap skipped — missing %s', f)
            return None

    df = pd.read_csv(gsf_csv)
    ctrl = list(pd.read_csv(flux1_csv, index_col=0, nrows=0).columns)
    test = list(pd.read_csv(flux2_csv, index_col=0, nrows=0).columns)
    ctrl_set, test_set = set(ctrl), set(test)
    sample_cols = [c for c in df.columns if c in ctrl_set or c in test_set]
    if not sample_cols:
        logging.warning('Pathway flux heatmap skipped — no sample flux columns found.')
        return None
    cond = {c: ('Control' if c in ctrl_set else 'Test') for c in sample_cols}

    # Enriched subsystems (pathway selection)
    if pathways:
        enriched = set(pathways)
    else:
        enriched = set()
        for fn in ('Up_regulated_pathways.csv', 'Down_regulated_pathways.csv'):
            p = os.path.join(output_dir, fn)
            if os.path.exists(p):
                pe = pd.read_csv(p)
                if 'Pathway' in pe.columns:
                    # Match the barplot (pathway_enrichment_plot.R), which keeps
                    # only pathways with Adjusted P-value < fdr_cutoff. Without this
                    # the heatmap legend would show CSV pathways that fail FDR and
                    # are therefore absent from the enrichment barplot.
                    if 'Adjusted P-value' in pe.columns:
                        pe = pe[pd.to_numeric(pe['Adjusted P-value'],
                                              errors='coerce') < fdr_cutoff]
                    enriched |= set(pe['Pathway'].dropna())
        if not enriched:
            logging.warning('Pathway flux heatmap skipped — no enriched pathways '
                            'pass FDR < %g.', fdr_cutoff)
            return None

    sub = df[df['Subsystem'].isin(enriched)].copy()
    if sub.empty:
        logging.warning('Pathway flux heatmap skipped — no rows in enriched pathways.')
        return None

    # Optional DE-flux-only restriction: keep only differentially-active reactions
    # (those in Differential_fluxes.csv), within the enriched pathways.
    if de_only:
        de_csv = os.path.join(output_dir, 'Differential_fluxes.csv')
        if not os.path.exists(de_csv):
            logging.warning('DE-only pathway heatmap skipped — missing %s', de_csv)
            return None
        de_rxns = set(pd.read_csv(de_csv, index_col=0).index.astype(str))
        sub = sub[sub['Reaction'].astype(str).isin(de_rxns)]
        if sub.empty:
            logging.warning('DE-only pathway heatmap skipped — no DE-flux reactions '
                            'in the enriched pathways.')
            return None

    # Keep the top-N most differentially-active features per pathway (by
    # |Effect size|, ties by ascending P). Caps figure size; keeps it labeled.
    if top_n_per_pathway:
        key = 'Reaction' if level == 'reaction' else 'Gene'
        score = pd.to_numeric(sub.get('Effect size'), errors='coerce').abs()
        if score.isna().all():
            score = pd.to_numeric(sub.get('log2 (condition2/condition1)'),
                                  errors='coerce').abs()
        pval = pd.to_numeric(sub.get('P'), errors='coerce')
        rank_df = sub[['Subsystem', key]].assign(
            _score=score.fillna(-np.inf).values, _p=pval.fillna(np.inf).values)
        # one score per (pathway, feature) = its strongest reaction
        feat = (rank_df.groupby(['Subsystem', key])
                       .agg({'_score': 'max', '_p': 'min'})
                       .reset_index()
                       .sort_values(['Subsystem', '_score', '_p'],
                                    ascending=[True, False, True], kind='mergesort'))
        feat['_rank'] = feat.groupby('Subsystem').cumcount() + 1
        keep = feat.loc[feat['_rank'] <= int(top_n_per_pathway), ['Subsystem', key]]
        sub = sub.merge(keep, on=['Subsystem', key], how='inner')
        if sub.empty:
            logging.warning('Pathway flux heatmap skipped — empty after top-N filter.')
            return None

    # HGNC symbol map (fallback to gene id)
    entrez2symbol = {}
    gm_file = _effective_gene_map(output_viz_dir)
    if os.path.exists(gm_file):
        gm = pd.read_csv(gm_file, dtype=str)
        if {'entrez', 'symbol'}.issubset(gm.columns):
            entrez2symbol = dict(zip(gm['entrez'], gm['symbol']))

    def _sym(g):
        s = entrez2symbol.get(str(g))
        return s if isinstance(s, str) and s and s.lower() != 'nan' else str(g)

    # ---- Build the feature × sample matrix + row annotations ------------------
    if level == 'reaction':
        sub_r = sub.drop_duplicates(subset='Reaction').set_index('Reaction')
        mat      = sub_r[sample_cols].apply(pd.to_numeric, errors='coerce')
        row_sub  = sub_r['Subsystem']
        row_stat = sub_r['Status']
        labels   = [str(r) for r in mat.index]
    else:  # gene
        g = sub.groupby('Gene')
        mat      = g[sample_cols].mean()                     # mean over the gene's reactions
        row_sub  = g['Subsystem'].agg(lambda s: s.mode().iloc[0])
        row_stat = g['Status'].agg(lambda s: s.mode().iloc[0])
        labels   = [_sym(gid) for gid in mat.index]

    mat = mat.fillna(0.0)
    label_of = dict(zip(mat.index, labels))

    # Z-score rows; drop zero-variance rows
    keep = mat.std(axis=1) > 0
    mat = mat.loc[keep]
    if mat.empty:
        logging.warning('Pathway flux heatmap skipped — all rows have zero variance.')
        return None
    row_sub, row_stat = row_sub.loc[mat.index], row_stat.loc[mat.index]
    Z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0).clip(-4, 4)

    # Columns ordered & split by Condition (Control then Test)
    col_order = ([c for c in sample_cols if cond[c] == 'Control']
                 + [c for c in sample_cols if cond[c] == 'Test'])
    Z = Z[col_order]

    # Display labels on the rows
    Z.index = [label_of[i] for i in Z.index]

    # ---- Annotation colors -----------------------------------------------------
    subs_u   = sorted(row_sub.unique())
    # colorblind-safe 10-color palette when few pathways; distinct hues (husl)
    # when many (>10), where keeping colors distinguishable matters most.
    _pal = 'colorblind' if len(subs_u) <= 10 else 'husl'
    sub_pal  = dict(zip(subs_u, sns.color_palette(_pal, max(len(subs_u), 1))))
    STAT_PAL = {'UP': '#D55E00', 'DOWN': '#0072B2', 'Reversed': '#666666',
                'Unchanged': '#CCCCCC'}
    COND_PAL = {'Control': '#009E73', 'Test': '#D55E00'}   # matches flux_dimred.R

    # Positional arrays (row_sub/row_stat are in Z row order)
    Zv          = Z.values
    row_sub_a   = list(row_sub.values)
    row_stat_a  = list(row_stat.values)
    row_labels  = list(Z.index)
    cols        = list(Z.columns)
    N, M        = Zv.shape

    # ---- Independent within-group clustering ----------------------------------
    # Rows clustered independently within each Direction (UP, DOWN, Reversed);
    # columns clustered independently within each Condition (Control, Test).
    from scipy.cluster.hierarchy import linkage, dendrogram
    from matplotlib.colors import to_rgb, Normalize
    from matplotlib.cm import ScalarMappable

    def _order(block):
        if block.shape[0] < 2:
            return list(range(block.shape[0])), None
        L = linkage(block, method='average', metric='euclidean')
        return dendrogram(L, no_plot=True)['leaves'], L

    row_order, row_groups = [], []          # row_groups: (status, count, linkage)
    for st in [s for s in ('UP', 'DOWN', 'Reversed') if s in set(row_stat_a)]:
        idx = [i for i in range(N) if row_stat_a[i] == st]
        leaves, L = _order(Zv[idx])
        row_groups.append((st, len(idx), L))
        row_order.extend(idx[j] for j in leaves)

    col_order_idx, col_groups = [], []      # col_groups: (condition, count, linkage)
    for cg in ('Control', 'Test'):
        cidx = [i for i in range(M) if cond[cols[i]] == cg]
        if not cidx:
            continue
        leaves, L = _order(Zv[:, cidx].T)
        col_groups.append((cg, len(cidx), L))
        col_order_idx.extend(cidx[j] for j in leaves)

    Zo         = Zv[np.ix_(row_order, col_order_idx)]
    rlab_o     = [row_labels[i] for i in row_order]
    clab_o     = [cols[i] for i in col_order_idx]
    rsub_o     = [row_sub_a[i] for i in row_order]
    rstat_o    = [row_stat_a[i] for i in row_order]
    ccond_o    = [cond[cols[i]] for i in col_order_idx]

    # ---- Figure layout (fixed fractional axes; legends + cbar all on the right)-
    set_publication_style()
    show_rlab = N <= 110
    fig_w = max(11.0, 0.13 * M + 7.0)
    fig_h = max(6.5, min(0.135 * N + 3.5, 46.0))
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')
    fig.set_constrained_layout(False)

    HX0, HW, HY0, HH = 0.165, 0.55, 0.07, 0.74           # heatmap rectangle
    RDx0, RDw = 0.030, 0.095                              # row dendrograms
    PAx, DAx, AW = 0.130, 0.145, 0.012                   # Pathway / Direction strips
    CAy, CAh = HY0 + HH + 0.004, 0.012                    # Condition strip
    CDy0, CDh = CAy + CAh + 0.004, 0.085                  # column dendrograms

    ax_hm = fig.add_axes([HX0, HY0, HW, HH])
    ax_hm.imshow(Zo, aspect='auto', cmap='RdBu_r', vmin=-4, vmax=4,
                 origin='upper', interpolation='nearest')
    ax_hm.set_xticks(np.arange(M)); ax_hm.set_yticks(np.arange(N) if show_rlab else [])
    ax_hm.set_xticklabels(clab_o, rotation=90, fontsize=8)
    if show_rlab:
        ax_hm.set_yticklabels(rlab_o, fontsize=8)
        ax_hm.yaxis.set_ticks_position('right')
    ax_hm.tick_params(length=0)
    for s in ax_hm.spines.values():
        s.set_visible(False)

    # Row dendrograms: one per Direction group, aligned to its row band (top→down)
    cum = 0
    for st, cnt, L in row_groups:
        top = HY0 + HH * (1 - cum / N)
        bot = HY0 + HH * (1 - (cum + cnt) / N)
        if L is not None:
            axd = fig.add_axes([RDx0, bot, RDw, top - bot]); axd.set_axis_off()
            dendrogram(L, ax=axd, orientation='left', no_labels=True,
                       color_threshold=0, above_threshold_color='#888888')
            axd.invert_yaxis()                            # leaf 0 → top (matches heatmap)
        cum += cnt

    # Column dendrograms: one per Condition group, aligned to its column band
    cum = 0
    for cg, cnt, L in col_groups:
        left  = HX0 + HW * (cum / M)
        right = HX0 + HW * ((cum + cnt) / M)
        if L is not None:
            axd = fig.add_axes([left, CDy0, right - left, CDh]); axd.set_axis_off()
            dendrogram(L, ax=axd, orientation='top', no_labels=True,
                       color_threshold=0, above_threshold_color='#888888')
        cum += cnt

    # Annotation strips (Pathway, Direction on the left; Condition on top)
    def _strip(x0, colors_rgb, vertical=True):
        ax = fig.add_axes([x0, HY0, AW, HH]) if vertical else \
             fig.add_axes([HX0, CAy, HW, CAh])
        arr = (np.array(colors_rgb).reshape(-1, 1, 3) if vertical
               else np.array(colors_rgb).reshape(1, -1, 3))
        ax.imshow(arr, aspect='auto', origin='upper'); ax.set_axis_off()
    _strip(PAx, [to_rgb(sub_pal[s]) for s in rsub_o])
    _strip(DAx, [to_rgb(STAT_PAL.get(s, '#CCCCCC')) for s in rstat_o])
    _strip(None, [to_rgb(COND_PAL[c]) for c in ccond_o], vertical=False)

    # Narrow flux z-score colorbar (right)
    cax = fig.add_axes([0.925, HY0, 0.010, HH * 0.30])
    cb = fig.colorbar(ScalarMappable(norm=Normalize(-4, 4), cmap='RdBu_r'),
                      cax=cax, ticks=[-4, -2, 0, 2, 4])
    cb.set_label('flux z-score', fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # Legends — all stacked on the right, with adaptive spacing so the (possibly
    # long) Pathway legend never overlaps the Direction/Condition legends.
    path_h = [Patch(color=sub_pal[s], label=textwrap.shorten(str(s), 30)) for s in subs_u]
    dir_h  = [Patch(color=STAT_PAL[k], label=k) for k in ('UP', 'DOWN', 'Reversed')
              if k in set(row_stat_a)]
    cond_h = [Patch(color=COND_PAL[k], label=k) for k in ('Control', 'Test')]
    entry = (7.5 * 1.7) / (72.0 * fig_h)                  # ~1 legend entry, fig-fraction
    top   = HY0 + HH
    dir_y  = top - entry * (len(path_h) + 2.0)
    cond_y = dir_y - entry * (len(dir_h) + 2.0)
    fig.legend(handles=path_h, title='Pathway', loc='upper left',
               bbox_to_anchor=(0.80, top), fontsize=7.5, title_fontsize=8.5, frameon=False)
    fig.legend(handles=dir_h, title='Direction', loc='upper left',
               bbox_to_anchor=(0.80, dir_y), fontsize=7.5, title_fontsize=8.5, frameon=False)
    fig.legend(handles=cond_h, title='Condition', loc='upper left',
               bbox_to_anchor=(0.80, cond_y), fontsize=7.5, title_fontsize=8.5, frameon=False)

    cohort = os.path.basename(os.path.normpath(output_dir))
    scope = 'DE-flux only' if de_only else 'all flux'
    fig.suptitle('%s — enriched-pathway flux (%s level, %s, z-scored)'
                 % (cohort, level, scope), x=0.45, y=0.965, fontsize=13,
                 fontweight='medium')

    suffix  = '_de' if de_only else ''
    out_png = os.path.join(output_viz_dir, 'pathway_flux_heatmap_%s%s.png' % (level, suffix))
    out_svg = os.path.join(output_viz_dir, 'pathway_flux_heatmap_%s%s.svg' % (level, suffix))
    os.makedirs(output_viz_dir, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_svg, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    logging.info('Saved: %s (+ .svg)  [%d %ss x %d samples, %d pathways, %s]',
                 out_png, N, level, M, len(subs_u), scope)
    return out_png


def _pathway_barplot_matplotlib(output_dir, output_viz_dir, fdr_cutoff=0.1, title=None):
    """Matplotlib backend for visualize_pathway_enrichment (use_r=False).

    Draws a single combined bar plot of -log10(FDR) for up- and down-regulated
    metabolic pathways, mirroring the R figure (same colors, ratio labels,
    publication style). Writes ``{output_viz_dir}/pathway_barplot.png`` (+ .svg)
    — the same output as the R backend.
    """
    import textwrap

    set_publication_style()
    os.makedirs(output_viz_dir, exist_ok=True)

    UP_COLOR, DOWN_COLOR = '#E41A1C', '#377EB8'

    def _load(direction):
        path = f'{output_dir}/{direction}_regulated_pathways.csv'
        if not os.path.exists(path):
            return None
        df = pd.read_csv(path, index_col=0)
        if df.empty:
            return None
        df = df.copy()
        df['Pathway']   = df.index
        df['Direction'] = direction.capitalize()   # 'Up' / 'Down'
        return df

    parts = [d for d in (_load('Up'), _load('Down')) if d is not None]
    if not parts:
        logging.warning('No pathway enrichment results found — skipping bar plot.')
        return

    combined = pd.concat(parts, ignore_index=True)
    combined = combined[combined['Adjusted P-value'] < fdr_cutoff]
    if combined.empty:
        logging.warning('No pathways pass FDR < %s — skipping bar plot.', fdr_cutoff)
        return

    combined['neg_log10_p'] = -np.log10(combined['Adjusted P-value'])
    overlap = combined['Overlap'].astype(str).str.split('/', expand=True)
    combined['k']       = overlap[0].astype(float)
    combined['n_total'] = overlap[1].astype(float)
    combined['enrich_r'] = (combined['k'] / combined['n_total']).round(2)
    combined['bar_label'] = (combined['enrich_r'].astype(str) + ' (' +
                             combined['k'].astype(int).astype(str) + '/' +
                             combined['n_total'].astype(int).astype(str) + ')')

    # Order: Up group then Down group, each ascending by significance (R: arrange(Direction, neg_log10_p))
    combined['Direction'] = pd.Categorical(combined['Direction'], categories=['Up', 'Down'], ordered=True)
    combined = combined.sort_values(['Direction', 'neg_log10_p']).reset_index(drop=True)

    n = len(combined)
    y_pos  = np.arange(n)
    colors = [UP_COLOR if d == 'Up' else DOWN_COLOR for d in combined['Direction']]
    fig_height = max(n * 0.45 + 2.8, 4)

    fig, ax = plt.subplots(figsize=(9, fig_height))
    ax.barh(y_pos, combined['neg_log10_p'], color=colors,
            edgecolor='black', linewidth=0.5, height=0.7)
    for yi, (val, lab) in enumerate(zip(combined['neg_log10_p'], combined['bar_label'])):
        ax.text(val / 2, yi, lab, ha='center', va='center', color='white', fontsize=10)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([textwrap.fill(p, 45) for p in combined['Pathway']])
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_xlim(0, combined['neg_log10_p'].max() * 1.08)
    ax.set_xlabel(r'$-\log_{10}$(FDR)')
    ax.set_ylabel('Metabolic pathway')
    ax.set_title(title or 'Enriched metabolic pathways', pad=30)

    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=UP_COLOR, edgecolor='black', label='Up'),
                       Patch(facecolor=DOWN_COLOR, edgecolor='black', label='Down')],
              title='Direction', loc='lower center', bbox_to_anchor=(0.5, 1.0),
              ncol=2, frameon=False)

    fig.savefig(f'{output_viz_dir}/pathway_barplot.png', format='png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_viz_dir}/pathway_barplot.svg', format='svg', bbox_inches='tight')
    plt.close(fig)
    logging.info('Saved: %s/pathway_barplot.png (+ .svg)', output_viz_dir)
    return
