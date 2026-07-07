import argparse
import copy
import logging
import os
import time
import warnings

import pandas as pd
from cobra.io import read_sbml_model
from tf_metabolism.metabolic_model import GPR_manipulation
from tf_metabolism.omics_integration import omics_score_calculation
from tf_metabolism.omics_integration import LAD


def parse_getpra_information(getpra_file):
    getrpa_sl_info = {}
    fp = open(getpra_file, 'r')
    fp.readline()
    for line in fp:
        sptlist = line.strip().split('\t')
        ucsc_id = sptlist[4].strip()
        compartment = sptlist[8].strip()

        if ucsc_id not in getrpa_sl_info:
            getrpa_sl_info[ucsc_id] = [compartment]
        else:
            getrpa_sl_info[ucsc_id].append(compartment)

    fp.close()

    return getrpa_sl_info


def calculate_reaction_score(cobra_model, expression_score, getrpa_sl_info, use_getpra_filtering=False):
    expression_socre_string = {}
    for gene_id in expression_score:
        expression_socre_string[str(gene_id)] = expression_score[gene_id]

    reaction_weights = {}
    for each_reaction in cobra_model.reactions:        
        temp_expression_socre_string = copy.deepcopy(expression_socre_string)

        compartments = [metabolite.compartment for metabolite in each_reaction.reactants + each_reaction.products]
        compartments = list(set(compartments))
        reaction_weights[each_reaction.id] = 0.0
        GPR_association = each_reaction.gene_reaction_rule
        genes = [gene.id for gene in each_reaction.genes]
        genes = list(set(genes) & set(temp_expression_socre_string.keys()))

        if use_getpra_filtering:
            for each_gene in genes:
                if each_gene in getrpa_sl_info:
                    getpra_compartments = getrpa_sl_info[each_gene]
                    if len(set(getpra_compartments) & set(compartments)) == 0:
                        temp_expression_socre_string[each_gene] = 0.0

        if len(genes) > 1:
            GPR_list = GPR_manipulation.convert_string_GPR_to_list_GPR(GPR_association)
            expression = omics_score_calculation.GPR_score_calculation(GPR_list, temp_expression_socre_string)
            reaction_weights[each_reaction.id] = expression
        elif len(genes) == 1:
            gene_id = genes[0].strip()
            if gene_id in temp_expression_socre_string:
                expression = temp_expression_socre_string[gene_id]
                reaction_weights[each_reaction.id] = expression

    return reaction_weights


def read_expressoin_data(omics_file, delimiter=','):
    expression_score = {}
    with open(omics_file, 'r') as fp:
        fp.readline()
        for line in fp:
            sptlist = line.strip().split(delimiter)
            isoform_id = sptlist[0].strip()
            expression_level = float(sptlist[1].strip())
            expression_score[isoform_id] = expression_level
    return expression_score

def calculate_flux(output_dir, omics_file, getpra_file, use_getpra, generic_cobra_model_file, context_specific_cobra_model_file, minimum_biomass=0.01, biomass_reaction_id='biomass_reaction'):
    # getpra_file is Recon2M.2-only (GeTPRA / BIGG). Human-GEM has none and passes None;
    # it is only consulted when use_getpra=True, so tolerate its absence. Default
    # biomass_reaction_id keeps the Recon2M.2 path byte-identical ('biomass_reaction').
    getrpa_sl_info = parse_getpra_information(getpra_file) if getpra_file else {}
    cobra_model = read_sbml_model(generic_cobra_model_file)
    cobra_isoforms = [isoform.id for isoform in cobra_model.genes]
    print('###', cobra_model)
    basename = os.path.basename(omics_file).split('.')[0].strip()
    gene_expression_df = pd.read_csv(omics_file, index_col=0)
    raw_expression_score = omics_score_calculation.calculate_rank_based_expression_score(gene_expression_df)
    expression_score = {k: max(v, 0.0) for k, v in raw_expression_score.items()}
    
    reaction_weight_info = {}
    if use_getpra == True:
        reaction_weights = calculate_reaction_score(cobra_model, expression_score, getrpa_sl_info, True)
    else:
        reaction_weights = calculate_reaction_score(cobra_model, expression_score, getrpa_sl_info, False)

    reaction_weight_info[basename] = reaction_weights

    df = pd.DataFrame.from_dict(reaction_weight_info)
    df.to_csv(output_dir + '/Reaction_weight_%s.csv' % (basename))
    
    context_specific_cobra_model = read_sbml_model(context_specific_cobra_model_file)
    
    obj = LAD.LAD()
    obj.load_cobra_model(context_specific_cobra_model)

    a,b,c = obj.run_FBA(new_objective=biomass_reaction_id)
    print('results : ', a, b)

    flux_constraints = {}
    flux_constraints[biomass_reaction_id] = [minimum_biomass, 1000.0]
    
    non_boundary_reactions = set(each_reaction.id for each_reaction in context_specific_cobra_model.reactions
                                 if not each_reaction.boundary)
    model_reactions = [each_reaction.id for each_reaction in context_specific_cobra_model.reactions]
    new_reaction_weight_info = {}
    for each_reaction in reaction_weights:
        if each_reaction in model_reactions:
            new_reaction_weight_info[each_reaction] = reaction_weights[each_reaction] # /10000.0
    solution_status, objective_value, predicted_flux = obj.run_LP_fitting(opt_flux=new_reaction_weight_info, flux_constraints=flux_constraints)
    print(solution_status, objective_value, predicted_flux)
    output_file = f'{output_dir}/Flux_prediction_{basename}.csv'
    print(predicted_flux)
    predicted_flux = {rxn: v for rxn, v in predicted_flux.items() if rxn in non_boundary_reactions}

    with open(output_file, 'w') as fp:
        for rxn, v in predicted_flux.items():
            fp.write(f'{rxn},{v}\n')

    return predicted_flux

def main():
    start = time.time()
    warnings.filterwarnings("ignore")

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_omics_data', required=True, help="Omics data")
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory")
    parser.add_argument('-getpra', '--getpra_file', required=True, help="GeTPRA file")
    parser.add_argument('-g_model', '--generic_cobra_model_file', required=True, help="Generic cobra model file")
    parser.add_argument('-c_model', '--context_specific_cobra_model_file', required=True,
                        help="Generic cobra model file")
    parser.add_argument('-use_getpra', '--use_getpra', default='yes', choices=['yes', 'no'], required=True,
                        help="Use GeTPRA")

    options = parser.parse_args()
    omics_file = options.input_omics_data
    output_dir = options.output_dir
    getpra_file = options.getpra_file
    use_getpra = options.use_getpra
    generic_cobra_model_file = options.generic_cobra_model_file
    context_specific_cobra_model_file = options.context_specific_cobra_model_file
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    calculate_flux(output_dir, omics_file, getpra_file, use_getpra, generic_cobra_model_file, context_specific_cobra_model_file)    

    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))


if __name__ == '__main__':
    main()
