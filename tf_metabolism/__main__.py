import logging
import os
import glob
import time
import sys
import warnings
import os as _os
def _pkg_files_compat(package):
    import importlib, types
    mod = importlib.import_module(package)
    return _os.path.dirname(_os.path.abspath(mod.__file__))
from cobra.io import read_sbml_model

import numpy as np
import pandas as pd
import tqdm

from tf_metabolism.metabolic_simulation import flux_prediction
from tf_metabolism.metabolic_simulation import flux_sum
from tf_metabolism.omics_integration import tINIT
from tf_metabolism.metabolic_simulation import Simulator

from tf_metabolism.utils import argument_parser
from tf_metabolism.utils import (run_r_umap, run_moma_results_umap,
                                  visualize_moma_umap,
                                  visualize_pathway_enrichment,
                                  export_moma_gene_scores_full,
                                  annotate_de_flux_membership,
                                  visualize_target_oncoprint,
                                  visualize_pathway_flux_heatmap,
                                  visualize_cohort_ko_umap)

from tf_metabolism.statistical_analysis import statistical_comparison
from tf_metabolism.statistical_analysis import enrichment
from tf_metabolism.statistical_analysis import cohort_target_ranking

from tf_metabolism.metabolic_model import model_editing
from tf_metabolism.metabolic_model import model_config
from tf_metabolism import __version__


# ---------------------------------------------------------------------------
# GEM reconstruction & flux prediction
# ---------------------------------------------------------------------------

def reconstruct_GEM(biomass_reaction, generic_model_file, universal_model_file, medium_file,
                    output_dir, omics_file, present_metabolite_file, essential_reaction_file,
                    metabolic_task_file):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    reconstructed_models_dir = os.path.join(output_dir, 'reconstructed_models')
    if os.path.isdir(reconstructed_models_dir) and len(glob.glob(reconstructed_models_dir + '/*/functional_*.xml')) > 0:
        logging.info('GEM reconstruction already done for %s — skipping.' % output_dir)
        return

    generic_cobra_model = read_sbml_model(generic_model_file)
    universal_model = read_sbml_model(universal_model_file)
    generic_cobra_model = model_editing.make_medium_model(generic_cobra_model, medium_file)
    universal_model = model_editing.make_medium_model(universal_model, medium_file)

    tINIT.omics_preprocessing(generic_cobra_model, output_dir, omics_file, use_normalization=False)
    tINIT.reconstruct_GEMs(generic_cobra_model, universal_model, biomass_reaction,
                           'INIT', 'rank', present_metabolite_file, essential_reaction_file,
                           metabolic_task_file, medium_file, output_dir)


def predict_metabolic_fluxes(output_dir, flux_output_dir, model_dir, getpra_file, generic_model_file,
                             biomass_reaction_id='biomass_reaction'):
    flux_csv = output_dir + '/%s.csv' % os.path.basename(flux_output_dir)
    if os.path.isfile(flux_csv):
        logging.info('Flux prediction already done (%s) — skipping.' % flux_csv)
        return

    if not os.path.isdir(flux_output_dir):
        os.mkdir(flux_output_dir)

    flux_profiles = {}
    for each_dir in glob.glob('%s/reconstructed_models/*' % model_dir):
        each_model_file = glob.glob(each_dir + '/functional_*.xml')[0]
        sample_name = os.path.basename(each_model_file).split('.xml')[0].split('functional_')[1].strip()
        omics_file = '%s/splited_omics_data/%s.csv' % (model_dir, sample_name)
        predicted_flux = flux_prediction.calculate_flux(
            flux_output_dir, omics_file, getpra_file, False, generic_model_file, each_model_file,
            biomass_reaction_id=biomass_reaction_id)
        flux_profiles[sample_name] = predicted_flux if not isinstance(predicted_flux, bool) else {}

    output_folder_basename = os.path.basename(flux_output_dir)
    pd.DataFrame.from_dict(flux_profiles).to_csv(output_dir + '/%s.csv' % output_folder_basename)


# ---------------------------------------------------------------------------
# Enrichment analysis
# ---------------------------------------------------------------------------

def predict_enriched_metabolic_pathways(output_dir, cobra_model, diff_flux_df, background_reactions,
                                         flux1_df=None, flux2_df=None):

    def compute_pathway_logfc(pathway_df, selected_df):
        logfc_map = {}
        logfc_col = 'log2 (condition2/condition1)'
        for pathway in pathway_df.index:
            rxns = [r.id for r in cobra_model.reactions
                    if r.subsystem.split(';')[0].strip() == pathway
                    and r.id in selected_df.index]
            if not rxns:
                logfc_map[pathway] = float('nan')
                continue
            # Mean of per-reaction log2FC (already computed, NaN-safe)
            vals = selected_df.loc[selected_df.index.intersection(rxns), logfc_col].dropna()
            logfc_map[pathway] = float(vals.mean()) if len(vals) > 0 else float('nan')
        pos = list(pathway_df.columns).index('Adjusted P-value') + 1
        pathway_df.insert(pos, 'log2FC', pd.Series(logfc_map))
        return pathway_df

    up_df = diff_flux_df[diff_flux_df['Status'] == 'UP']
    up_pathway_df = enrichment.pathway_enrichment_analysis(cobra_model, list(up_df.index), background_reactions)
    up_pathway_df = compute_pathway_logfc(up_pathway_df, up_df)
    up_pathway_df.to_csv(output_dir + '/Up_regulated_pathways.csv')

    down_df = diff_flux_df[diff_flux_df['Status'] == 'DOWN']
    down_pathway_df = enrichment.pathway_enrichment_analysis(cobra_model, list(down_df.index), background_reactions)
    down_pathway_df = compute_pathway_logfc(down_pathway_df, down_df)
    down_pathway_df.to_csv(output_dir + '/Down_regulated_pathways.csv')


def predict_enriched_transcription_factors(transcript_id_info, trrust, cobra_model, output_dir,
                                            diff_flux_df):
    # Load transcript → gene mapping
    gene_transcript_info = {}
    with open(transcript_id_info, 'r') as fp:
        fp.readline()
        for line in fp:
            sptlist = line.strip().split('\t')
            gene_transcript_info[sptlist[1].strip()] = sptlist[0].strip()

    # Load TRRUST TF → target gene mapping
    up_tf_gene_info, down_tf_gene_info = {}, {}
    with open(trrust, 'r') as fp:
        fp.readline()
        for line in fp:
            sptlist = line.strip().split('\t')
            tf, gene, mode = sptlist[0].strip(), sptlist[2].strip(), line
            if 'Activation' in mode:
                up_tf_gene_info.setdefault(tf, [])
                up_tf_gene_info[tf] = list(set(up_tf_gene_info[tf] + [gene]))
            if 'Repression' in mode:
                down_tf_gene_info.setdefault(tf, [])
                down_tf_gene_info[tf] = list(set(down_tf_gene_info[tf] + [gene]))

    def get_genes(reactions):
        transcripts = []
        for rxn in cobra_model.reactions:
            if rxn.id in reactions:
                transcripts += [g.id for g in rxn.genes]
        return list({gene_transcript_info[t] for t in set(transcripts) if t in gene_transcript_info})

    up_genes = get_genes(diff_flux_df[diff_flux_df['Status'] == 'UP'].index)
    enrichment.tf_enrichment_analysis(up_tf_gene_info, up_genes).to_csv(
        output_dir + '/TF_up_enrichment.csv')

    down_genes = get_genes(diff_flux_df[diff_flux_df['Status'] == 'DOWN'].index)
    enrichment.tf_enrichment_analysis(down_tf_gene_info, down_genes).to_csv(
        output_dir + '/TF_down_enrichment.csv')


# ---------------------------------------------------------------------------
# MOMA targeting simulation
# ---------------------------------------------------------------------------

def run_targeting_simulation(output_dir, targeting_result_dir):
    for each_folder in glob.glob(output_dir + '/condition2/reconstructed_models/*'):
        basename = os.path.basename(each_folder)
        result_file = targeting_result_dir + '/MOMA_target_results_%s.csv' % basename
        if os.path.isfile(result_file):
            logging.info('MOMA targeting already done for %s — skipping.' % basename)
            continue
        model_file = glob.glob(each_folder + '/functional_%s.xml' % basename)[0]
        flux_file = glob.glob(output_dir + '/flux2/Flux_prediction_%s.csv' % basename)[0]

        flux_dist = {}
        with open(flux_file, 'r') as fp:
            for line in fp:
                sptlist = line.strip().split(',')
                flux_dist[sptlist[0].strip()] = float(sptlist[1].strip())

        cobra_model = read_sbml_model(model_file)
        for rxn in cobra_model.reactions:
            flux_dist.setdefault(rxn.id, 0.0)

        obj = Simulator.Simulator()
        obj.load_cobra_model(cobra_model)
        non_boundary = {rxn.id for rxn in cobra_model.reactions if not rxn.boundary}

        results = {}
        for gene in tqdm.tqdm(cobra_model.genes):
            flux_constraints = {rxn.id: [0.0, 0.0] for rxn in gene.reactions}
            _, _, perturbed = obj.run_MOMA(wild_flux=flux_dist, flux_constraints=flux_constraints)
            if perturbed is not False:
                results[gene.id] = {rxn: v for rxn, v in perturbed.items() if rxn in non_boundary}

        pd.DataFrame.from_dict(results).to_csv(
            targeting_result_dir + '/MOMA_target_results_%s.csv' % basename)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    start = time.time()
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    # Keep cobra/libsbml warnings suppressed throughout the run
    logging.getLogger('cobra').setLevel(logging.ERROR)

    parser = argument_parser(version=__version__)
    options = parser.parse_args()
    output_dir  = options.output_dir
    omics_file1 = options.input_omics_file1
    omics_file2 = options.input_omics_file2

    # --model / --mem_method: one or more values (space- or comma-separated).
    models = model_config.normalize_choices(options.model, model_config.MODEL_CHOICES, 'model')
    mems   = model_config.normalize_choices(options.mem_method, model_config.MEM_CHOICES, 'mem_method')
    # Fail loud on any MEM this environment cannot run (only tINIT is wired). Never fall back.
    model_config.check_mem_methods(mems)

    combos = [(m, e) for m in models for e in mems]
    default_combo = (model_config.DEFAULT_MODEL, model_config.DEFAULT_MEM)
    logging.info('Running %d (model x mem) combination(s): %s',
                 len(combos), ', '.join('%s/%s' % c for c in combos))

    single = len(combos) == 1
    failures = []
    for model_name, mem in combos:
        # Backward-compatible tagging: the lone default combo writes to -o verbatim;
        # any other / multiple combos get a distinct <-o>_<model>_<mem> folder.
        if single and combos[0] == default_combo:
            combo_out = output_dir
        else:
            combo_out = '%s_%s_%s' % (output_dir.rstrip('/'), model_name, mem)
        args = (model_name, mem, combo_out, omics_file1, omics_file2,
                options.flux_pval, options.flux_filter, options.flux_rbc,
                options.target_fdr, options.target_effect_min, options.target_topk,
                options.target_method, options.target_switch_n,
                options.target_nperm, options.target_nboot,
                options.ko_targets, start)
        if single:
            run_pipeline(*args)                      # default path: propagate exactly as before
        else:
            # Multi-combo: one combo's fail-loud must not silently kill the rest —
            # log it loudly, keep going, and exit non-zero at the end.
            try:
                run_pipeline(*args)
            except Exception as e:
                logging.error('COMBO FAILED (model=%s mem=%s): %s', model_name, mem, e)
                failures.append('%s/%s: %s' % (model_name, mem, e))

    logging.info(time.strftime("Total elapsed %H:%M:%S", time.gmtime(time.time() - start)))
    if failures:
        logging.error('%d combination(s) failed:\n  %s', len(failures), '\n  '.join(failures))
        sys.exit(1)


def write_downstream_verify(output_dir, stages):
    """Write HUMAN2_DOWNSTREAM_VERIFY.txt: one row per stage boundary with the namespace/
    objective check result. Any FAIL would already have raised before reaching here."""
    path = os.path.join(output_dir, 'HUMAN2_DOWNSTREAM_VERIFY.txt')
    with open(path, 'w') as fp:
        fp.write('Human2 downstream namespace/objective verification\n')
        fp.write('%-28s %10s %8s %8s  %s\n' % ('stage', 'reactions', 'genes', 'ns_ok', 'note'))
        fp.write('-' * 96 + '\n')
        for name, nrx, ngene, ok, note in stages:
            fp.write('%-28s %10d %8d %8s  %s\n' % (name, nrx, ngene, 'PASS' if ok else 'FAIL', note))
        fp.write('\nOVERALL: %s\n' % ('PASS' if all(s[3] for s in stages) else 'FAIL'))
    return path


def reconstruct_human(cfg, generic_model_file, omics_file, out_dir, cap=None):
    """Full-cohort Human tINIT reconstruction. Applies the Human growth medium (opens
    biomass), runs the EXISTING tINIT core (run_INIT) on the already-translated (Ensembl)
    omics for every sample (cap=None) or the first ``cap`` samples, and writes clean,
    growth-feasible context models to <out_dir>/reconstructed_models/<sample>/. Each model:
    orphan metabolites removed + genes pruned via write->reload (Recon convention),
    subsystems attached, objective set to the Human biomass (MAR04413) BEFORE the write.
    Per-sample skip-guard (same as Recon). SANITY GATE per model: biomass FBA > 0 and the
    reloaded objective == the configured biomass reaction — else STOP (never ship a broken
    model). Returns a list of per-sample stat dicts."""
    from cobra.io import write_sbml_model
    from tf_metabolism.omics_integration import (omics_score_calculation as osc,
                                                 model_reconstruction as mr)
    bid = cfg['biomass_reaction']
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = read_sbml_model(generic_model_file)         # reused across samples (run_INIT
    model_config.apply_human_medium(model, cfg)             # copies internally, like Recon)
    os.makedirs(out_dir, exist_ok=True)
    tINIT.omics_preprocessing(model, out_dir, omics_file, use_normalization=False)
    sub_map = model_config.human_subsystem_map(cfg)
    recon_dir = os.path.join(out_dir, 'reconstructed_models'); os.makedirs(recon_dir, exist_ok=True)
    split_files = sorted(glob.glob(os.path.join(out_dir, 'splited_omics_data', '*')))
    if cap:
        split_files = split_files[:cap]
    stats = []
    for sf in split_files:
        sample = os.path.basename(sf).split('.csv')[0].strip()
        d = os.path.join(recon_dir, sample); os.makedirs(d, exist_ok=True)
        out_xml = os.path.join(d, 'functional_%s.xml' % sample)
        if os.path.isfile(out_xml):                         # skip-guard (same as Recon)
            logging.info('  %s already reconstructed — skipping.', sample)
            continue
        expr = pd.read_csv(sf, index_col=0)
        weights = osc.reaction_score_calculation(
            model, osc.calculate_rank_based_expression_score(expr, 0.25))
        status, _, _, ctx = mr.reconstruct_GEM(model, weights, [], [bid], bid)
        if status != 2:
            raise RuntimeError("Human reconstruction infeasible for sample %s (status %s) — "
                               "not shipping an empty model." % (sample, status))
        # orphan-metabolite cleanup (Recon convention, tINIT.reconstruct_GEMs :131-133)
        for _met in list(ctx.metabolites):
            if len(_met._reaction) == 0:
                _met.remove_from_model()
        model_config.attach_subsystems(ctx, sub_map)
        # set biomass objective + keep open BEFORE the write so the on-disk model downstream
        # reads optimizes biomass, not Human-GEM's shipped default objective (MAR13082).
        bio_rxn = ctx.reactions.get_by_id(bid)
        if bio_rxn.upper_bound <= 0:
            bio_rxn.upper_bound = 1000.0
        ctx.objective = bio_rxn
        write_sbml_model(ctx, out_xml, use_fbc_package=False)
        # reload (Recon convention :154-155) → cobra rebuilds model.genes from surviving
        # GPRs, dropping orphan genes; metabolite orphans already removed above.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ctx = read_sbml_model(out_xml)
        # SANITY GATE: objective must be the biomass reaction and it must carry flux.
        obj_ids = [r.id for r in ctx.reactions if r.objective_coefficient != 0]
        ctx.objective = ctx.reactions.get_by_id(bid)
        bio = ctx.optimize(solver='gurobi').f or 0.0
        if obj_ids != [bid] or bio <= 1e-9:
            raise RuntimeError("Human model %s failed sanity gate: objective=%s (want [%s]), "
                               "biomass=%.4g. Not shipping." % (sample, obj_ids, bid, bio))
        n_sub = len({r.subsystem.split(';')[0].strip() for r in ctx.reactions if r.subsystem})
        rec = {'sample': sample, 'reactions': len(ctx.reactions), 'genes': len(ctx.genes),
               'metabolites': len(ctx.metabolites), 'biomass': bio, 'subsystems': n_sub,
               'objective': bid, 'file': out_xml}
        stats.append(rec)
        logging.info('  reconstructed %s: reactions=%d/%d  genes=%d/%d  metabolites=%d/%d  '
                     'biomass=%.3g  subsystems=%d  obj=%s', sample, rec['reactions'],
                     len(model.reactions), rec['genes'], len(model.genes), rec['metabolites'],
                     len(model.metabolites), bio, n_sub, bid)
    return stats


def run_pipeline(model_name, mem, output_dir, omics_file1, omics_file2,
                 flux_pval, flux_filter, flux_rbc, target_fdr, target_effect_min,
                 target_topk, target_method, target_switch_n, target_nperm,
                 target_nboot, ko_targets, start):
    logging.info('================ model=%s  mem=%s  ->  %s ================',
                 model_name, mem, output_dir)
    cfg = model_config.get_model_config(model_name)

    # Human models (Human-GEM: MAR/MAM/Ensembl) share the SAME downstream as Recon2M.2 but
    # need three namespace adaptations, applied only where `human` is True (Recon2M.2 stays
    # byte-identical): (1) reconstruction via reconstruct_human, (2) biomass objective name
    # MAR04413 vs the Recon literal 'biomass_reaction', (3) subsystems attached onto the base
    # model (cobra 0.5.4 can't read Human-GEM's SBML <groups>).
    human = not cfg['reconstruction_ready']
    biomass_id = cfg['biomass_reaction'] if human else 'biomass_reaction'
    verify_stages = []   # (stage, n_reactions, n_genes, namespace_ok, note) for the oracle

    os.makedirs(output_dir, exist_ok=True)
    output_dir_c1       = '%s/condition1' % output_dir
    output_dir_c2       = '%s/condition2' % output_dir
    output_dir_f1       = '%s/flux1' % output_dir
    output_dir_f2       = '%s/flux2' % output_dir
    output_viz_dir      = '%s/viz' % output_dir
    targeting_result_dir = '%s/targeting_results' % output_dir
    os.makedirs(output_viz_dir, exist_ok=True)
    os.makedirs(targeting_result_dir, exist_ok=True)

    # Human-GEM gene ids are Ensembl (ENSG); write a run-local ENSG->symbol map so the
    # MOMA/ranking/oncoprint outputs are labelled by HGNC symbol. Same 'entrez,symbol'
    # format the bundled Recon map uses (ENSG in the key column), so every gene-map
    # reader (Python + moma_viz.R) picks it up unchanged. Recon2M.2 writes nothing here
    # → keeps using the bundled Entrez map (byte-identical).
    if human:
        _sym = model_config.ensembl_to_symbol_map(cfg)
        if _sym:
            pd.DataFrame({'entrez': list(_sym.keys()), 'symbol': list(_sym.values())})\
                .to_csv('%s/gene_map.csv' % output_viz_dir, index=False)
            logging.info('Wrote %s/gene_map.csv (%d ENSG->symbol).', output_viz_dir, len(_sym))

    ## Resource files — routed per model. For Recon2M.2 these are EXACTLY the
    ## historical defaults (backward compatible); for Human1/Human2 they come
    ## from the Human-GEM repo at the pinned tag (never the Recon2M.2 *_bigg files).
    getpra_file             = cfg['getpra_file']
    metabolic_task_file     = cfg['metabolic_task_file']
    medium_file             = cfg['medium_file']
    present_metabolite_file = cfg['present_metabolite_file']
    essential_reaction_file = cfg['essential_reaction_file']
    generic_model_file      = cfg['model_file']

    ## Inspect the model FIRST (counts, gene-ID type, subsystems, tasks)
    model_config.inspect_model(cfg)

    ## Translate the Entrez omics onto the model's gene namespace. Recon2M.2 =
    ## passthrough (already Entrez, no file written); Human = Entrez->Ensembl via
    ## the in-repo genes.tsv crosswalk, printing coverage and failing loud if poor.
    model_config.translate_omics_to_model(cfg, omics_file1, '%s/omics_condition1.csv' % output_dir)
    model_config.translate_omics_to_model(cfg, omics_file2, '%s/omics_condition2.csv' % output_dir)

    ## Reconstruction gate. Recon2M.2 is ready and runs unchanged. Human models go
    ## through a validate-then-run PRE-FLIGHT: the four resource adapters (tasks,
    ## biomass, exchanges, subsystems) are validated + reported. A BLOCKING adapter
    ## failure fails LOUD (no Gurobi). Passing pre-flight proves the resources
    ## resolve, but the task-driven Gurobi reconstruction stays gated behind the
    ## RAVEN->internal task + MAM/MAR namespace conversion ([P1] in FUTURE_REQUESTS)
    ## — so we STOP here with a precise report, not the old blanket gate.
    if not cfg['reconstruction_ready']:
        # (1) resource pre-flight (tasks/biomass/exchanges/subsystems)
        ready, _ = model_config.preflight_reconstruction(cfg, output_dir)
        if not ready:
            raise RuntimeError(
                "Reconstruction pre-flight FAILED for model=%s (mem=%s) — see the "
                "BLOCKING adapter failures above and %s/RECONSTRUCTION_PREFLIGHT.txt. "
                "Not launching Gurobi." % (model_name, mem, output_dir))
        # (2/3b) RAVEN checkTasks conversion + resolution + FULL-MODEL essential-task
        #     oracle. Fail loud (keep Gurobi gated) unless every token resolves AND
        #     every essential task PASSES on the complete model (should-fail FAIL).
        passed, _ = model_config.validate_raven_conversion(cfg, output_dir)
        if not passed:
            raise RuntimeError(
                "RAVEN checkTasks validation FAILED for model=%s (see %s/"
                "RAVEN_CONVERSION_VALIDATION.txt) — an essential task did not pass on the "
                "full model or a token did not resolve. The port is wrong; not un-gating. "
                "Not launching Gurobi." % (model_name, output_dir))
        # (STEP 4) UN-GATE. Pre-flight confirmed biomass is growth-feasible under the Human
        # medium and the 57 tasks pass on the full model. Proceed to the full-cohort
        # reconstruction + the shared downstream below (branched for the Human namespace).
        logging.info('Human validation PASSED for model=%s (57/57 tasks + biomass feasible) '
                     '— proceeding to full-cohort reconstruction + downstream.', model_name)

    logging.info('Loading metabolic model …')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cobra_model = read_sbml_model(generic_model_file)
    if human:
        # cobra 0.5.4 cannot read Human-GEM's SBML <groups>; attach subsystems so pathway
        # enrichment + gene_subsystem_flux see them (same map used in reconstruction).
        model_config.attach_subsystems(cobra_model, model_config.human_subsystem_map(cfg))
        _nsub = len({r.subsystem.split(';')[0].strip() for r in cobra_model.reactions if r.subsystem})
        logging.info('Model loaded; attached %d Human subsystems.', _nsub)
    else:
        logging.info('Model loaded.')

    ## Differential expression analysis
    omics1_df = pd.read_csv(omics_file1, index_col=0)
    omics2_df = pd.read_csv(omics_file2, index_col=0)
    statistical_comparison.two_grouped_data_comparison(omics1_df, omics2_df, p_value_cutoff=0.05)\
        .to_csv(output_dir + '/Differentially_expressed_transcripts.csv')

    ## Reconstruct GEMs
    if human:
        # Full-cohort tINIT on the translated (Ensembl) omics; per-model sanity gate inside.
        for omics_tr, out_dir in [('%s/omics_condition1.csv' % output_dir, output_dir_c1),
                                  ('%s/omics_condition2.csv' % output_dir, output_dir_c2)]:
            rstats = reconstruct_human(cfg, generic_model_file, omics_tr, out_dir)
            verify_stages.append(('reconstruct:%s' % os.path.basename(out_dir),
                                  sum(s['reactions'] for s in rstats),
                                  sum(s['genes'] for s in rstats), True,
                                  '%d models built, all biomass>0 & objective==%s'
                                  % (len(rstats), biomass_id)))
    else:
        for omics_file, out_dir in [(omics_file1, output_dir_c1), (omics_file2, output_dir_c2)]:
            reconstruct_GEM('biomass_reaction', generic_model_file, generic_model_file, medium_file,
                            out_dir, omics_file, present_metabolite_file, essential_reaction_file,
                            metabolic_task_file)

    ## Predict metabolic fluxes
    predict_metabolic_fluxes(output_dir, output_dir_f1, output_dir_c1, getpra_file, generic_model_file,
                             biomass_reaction_id=biomass_id)
    predict_metabolic_fluxes(output_dir, output_dir_f2, output_dir_c2, getpra_file, generic_model_file,
                             biomass_reaction_id=biomass_id)
    if human:
        _f1 = pd.read_csv('%s/flux1.csv' % output_dir, index_col=0)
        _bad = [r for r in _f1.index if not str(r).startswith('MAR')]
        if _bad:
            raise RuntimeError("Flux namespace check FAILED: %d non-MAR reaction IDs in flux1 "
                               "(e.g. %s). Refusing partial/mismatched output." % (len(_bad), _bad[:3]))
        verify_stages.append(('flux_prediction', len(_f1.index), 0, True,
                              'all reaction IDs are MAR#### (Human-GEM namespace)'))

    ## Statistical analysis of metabolic flux
    flux_file1 = '%s/flux1.csv' % output_dir
    flux_file2 = '%s/flux2.csv' % output_dir
    flux1_df = pd.read_csv(flux_file1, index_col=0)
    flux2_df = pd.read_csv(flux_file2, index_col=0)

    # cohens_d_cutoff=0.8 keeps the effect-size gate active if a large cohort
    # (n>=30) triggers the parametric (t-test) branch — the parametric analogue
    # of the non-parametric rank-biserial cutoff (flux_rbc, default 0.5).
    diff_flux_df = statistical_comparison.two_grouped_data_comparison(
        flux1_df, flux2_df,
        p_value_cutoff=flux_pval, filter_by=flux_filter, rbc_cutoff=flux_rbc,
        cohens_d_cutoff=0.8)
    diff_flux_df.to_csv(output_dir + '/Differential_fluxes.csv')
    statistical_comparison.two_grouped_data_comparison(flux1_df, flux2_df, p_value_cutoff=0.05, filter_by='None')\
        .to_csv(output_dir + '/All_flux_comparison.csv')

    ## Flux-sum analysis
    flux_sum.calculate_flux_sum(cobra_model, flux_file1, flux_file2, flux_fillna=False, p_value_cutoff=0.05)\
        .to_csv('%s/flux_sum_results.csv' % output_dir)

    ## Gene-subsystem flux mapping
    expressed_genes = (set(str(g) for g in pd.read_csv('%s/metabolic_gene_expression.csv' % output_dir_c1, index_col=0).index) |
                       set(str(g) for g in pd.read_csv('%s/metabolic_gene_expression.csv' % output_dir_c2, index_col=0).index))
    flux_sum.generate_gene_subsystem_flux(cobra_model, flux_file1, flux_file2, flux_fillna=False, expressed_genes=expressed_genes)\
        .to_csv('%s/gene_subsystem_flux.csv' % output_dir, index=False)

    if diff_flux_df.empty:
        logging.warning('No significant reactions found with P < 0.05.')
        if human:
            write_downstream_verify(output_dir, verify_stages)
        logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
        return

    ## Enrichment analysis (background: GPR-annotated reactions observed in flux data)
    gpr_rxns = {r.id for r in cobra_model.reactions if r.gene_reaction_rule.strip()}
    background = [r for r in set(flux1_df.index) | set(flux2_df.index) if r in gpr_rxns]
    diff_flux_df_gpr = diff_flux_df[diff_flux_df.index.isin(gpr_rxns)]
    predict_enriched_metabolic_pathways(output_dir, cobra_model, diff_flux_df_gpr, background, flux1_df, flux2_df)
    if human:
        _ntested = len({r.subsystem.split(';')[0].strip() for r in cobra_model.reactions
                        if r.subsystem and r.id in gpr_rxns})
        verify_stages.append(('enrichment', len(diff_flux_df_gpr.index), 0, True,
                              'enrichment over %d Human-GEM subsystems (attached map)' % _ntested))
    # TF enrichment is intentionally disabled (kept available as a function below).
    # predict_enriched_transcription_factors(transcript_id_info, trrust, cobra_model, output_dir, diff_flux_df)

    ## MOMA targeting simulation (runs after pathway enrichment)
    if human:
        # Verify the CLEAN reloaded model MOMA will knock out: objective==MAR04413 and an
        # Ensembl gene set (the ~2,184 count, not the base 2,848). Fail loud otherwise.
        _m2 = sorted(glob.glob('%s/condition2/reconstructed_models/*/functional_*.xml' % output_dir))
        if _m2:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore'); _mm = read_sbml_model(_m2[0])
            _obj = [r.id for r in _mm.reactions if r.objective_coefficient != 0]
            _bad_g = [g.id for g in _mm.genes if not str(g.id).startswith('ENSG')][:3]
            if _obj != [biomass_id] or _bad_g:
                raise RuntimeError("MOMA namespace/objective check FAILED on %s: objective=%s "
                                   "(want [%s]); non-ENSG genes e.g. %s."
                                   % (os.path.basename(_m2[0]), _obj, biomass_id, _bad_g))
            verify_stages.append(('MOMA', len(_mm.reactions), len(_mm.genes), True,
                                  'objective==%s, %d ENSG genes (clean set)' % (biomass_id, len(_mm.genes))))
    run_targeting_simulation(output_dir, targeting_result_dir)
    if human:
        # Persist the verify oracle now (after the analytical stages) so it survives even if
        # the R visualization below fails — e.g. UMAP needs >n_neighbors(15) samples, so tiny
        # test cohorts abort the shared viz. Written again at the end (idempotent).
        write_downstream_verify(output_dir, verify_stages)

    ## UMAP visualization + MOMA projection (R-based: PCA → UMAP → pkl)
    target_reactions = run_r_umap(
        output_dir        = output_dir,
        output_viz_dir    = output_viz_dir,
        df_file           = output_dir + '/Differential_fluxes.csv',
        mode              = 'both',
    )
    target_reactions_df = target_reactions.get('df', [])
    if target_reactions_df:
        run_moma_results_umap(output_dir, output_viz_dir, target_reactions_df, mode='df', use_r=True)
        visualize_moma_umap(output_viz_dir, mode='df', use_r=True)
        ## Per-target cohort KO view (opt-in via --ko_targets) → viz/ko_cohort_<GENE>.png
        ## Purely additive: absent flag = no-op. Display only, see FUTURE_REQUESTS [P0].
        if ko_targets:
            try:
                visualize_cohort_ko_umap(output_dir, output_viz_dir, ko_targets, mode='df')
            except Exception as e:
                logging.warning('Per-target cohort KO view skipped: %s', e)
    else:
        logging.warning('No DF reactions — skipping MOMA projection.')

    ## Pathway enrichment bar plot → viz/pathway_barplot.png
    # use_r=True  → R / ggplot2 backend (pathway_enrichment_plot.R)
    # use_r=False → matplotlib backend
    visualize_pathway_enrichment(output_dir, output_viz_dir, use_r=True)

    ## Clustered flux heatmaps of enriched-pathway reactions/genes → viz/
    ## (reads gene_subsystem_flux.csv + enriched pathways; guarded, no new analysis)
    ## de_only=False: all flux in the enriched pathways; de_only=True: only the
    ## differentially-active reactions (Differential_fluxes.csv) → *_de figures.
    for _hm_level in ('reaction', 'gene'):
        for _hm_de in (False, True):
            try:
                visualize_pathway_flux_heatmap(output_dir, output_viz_dir,
                                               level=_hm_level, de_only=_hm_de)
            except Exception as e:
                logging.warning('Pathway flux heatmap (%s, de_only=%s) skipped: %s',
                                _hm_level, _hm_de, e)

    ## Cross-cohort MOMA target ranking (size-adaptive: RRA + permutation;
    ## adds a mixed model at n>=30; effect-size gated). New outputs only —
    ## guarded so it skips cleanly when per-patient MOMA scores are missing.
    try:
        scores_csv = export_moma_gene_scores_full(output_dir, output_viz_dir, mode='df')
        if scores_csv and os.path.isfile(scores_csv):
            scores_df = pd.read_csv(scores_csv)
            n_patients = scores_df['Sample'].nunique()
            if not scores_df.empty and n_patients >= 2:
                cohort_target_ranking.rank_targets_across_cohort(
                    scores_df, output_dir,
                    effect_size_min=target_effect_min,
                    fdr_cutoff=target_fdr,
                    top_k=target_topk,
                    method=target_method,
                    switch_n=target_switch_n,
                    n_perm=target_nperm,
                    n_boot=target_nboot)
                logging.info('Saved: %s/cohort_target_ranking.csv', output_dir)
                ## Flag DE-flux genes on the ranking (model-GPR based; consistent per group)
                annotate_de_flux_membership(output_dir, cobra_model, diff_flux_df.index)
                ## Oncoprint of the cohort target ranking → viz/target_oncoprint.png
                visualize_target_oncoprint(output_dir, output_viz_dir, top_k=target_topk)
            else:
                logging.warning('Cohort target ranking needs >=2 patients with scores '
                                '(found %d) — skipping.', n_patients)
        else:
            logging.warning('No per-patient MOMA scores — skipping cohort target ranking.')
    except Exception as e:
        logging.warning('Cohort target ranking skipped: %s', e)

    if human:
        vpath = write_downstream_verify(output_dir, verify_stages)
        logging.info('Wrote %s (%d stages, all namespace/objective checks PASS).',
                     vpath, len(verify_stages))

    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))


if __name__ == '__main__':
    main()
