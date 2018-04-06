import cobra
from cobra.flux_analysis import flux_variability_analysis
from sympy import Symbol
from cobrame.solve.algorithms import solve_at_growth_rate, binary_search
from cobrame import mu


def make_auxotrophs(model_cons, ko, rename_dict):
    model = model_cons.copy()
    if len(ko) > 0 and ko[0] in model_cons.reactions:
        for rxn in ko:
            model_cons.reactions.get_by_id(rxn).knock_out()
    else:
        ko_bnums = [rename_dict[i] for i in ko]
        cobra.manipulation.remove_genes(model, ko_bnums)
    model.id = str(ko)
    model.repair()
    return model


def _add_explicit_bound_constraint(model):
    for r in model.reactions:
        if r.id + '__UB_constraint' not in model.metabolites:
            met = cobra.Metabolite(r.id + '__UB_constraint')
            model.add_metabolites(met)
            met._bound = r.upper_bound
            met._constraint_sense = 'L'
            r.add_metabolites({met: 1.})
            r.upper_bound = 1000.

        if r.id + '__LB_constraint' not in model.metabolites:
            # if r.lower_bound != 0:
            met = cobra.Metabolite(r.id + '__LB_constraint')
            model.add_metabolites(met)
            met._bound = -r.lower_bound
            met._constraint_sense = 'L'
            r.add_metabolites({met: -1.})

        r.lower_bound = -1000
        r.upper_bound = 1000.


def _add_constraints_to_abundance_strain_variable(model, abundance_rxn):
    constraint_dict = {}
    for met in model.metabolites.query('_constraint'):
        coeff = -met._bound if '__UB_' in met.id else -met._bound
        constraint_dict[met] = coeff
        met._bound = 0
    abundance_rxn.add_metabolites(constraint_dict)

    biomass_dict = {}
    for biomass_rxn in model.objective.keys():
        for met, value in biomass_rxn.metabolites.items():
            if 'constraint' not in met.id:
                biomass_dict[met] = value * mu
        biomass_rxn.objective_coefficient = 0
    abundance_rxn.add_metabolites(biomass_dict)


def add_model_to_community(com_model, cons_model, i):

    model = cons_model.copy()

    print(i, model.id, model.optimize().f)
    abundance = cobra.Reaction('abundance')
    abundance.upper_bound = 1000.
    model.add_reaction(abundance)


    print(i, model.id, model.optimize().f)

    # add metabolites to new model with appended ID for strain compartment
    for met in model.metabolites:
        met.id += '_strain_%i' % i
    model.repair()

    for rxn in model.reactions:

        new_rxn = rxn.copy()
        new_rxn.id += '_strain_%i' % i

        com_model.add_reaction(new_rxn)

        if rxn.id.startswith('EX_'):
            coeff = 1 if 'reverse' not in rxn.id else -1
            shared_met = \
            rxn.id.replace('EX_', '').replace('_reverse', '').split('_strain')[
                0] + '_shared'
            new_rxn.add_metabolites({shared_met: coeff})
            #new_rxn.upper_bound = 1000

    _add_explicit_bound_constraint(model)

    abundance = com_model.reactions.get_by_id('abundance_strain_%i' % i)
    _add_constraints_to_abundance_strain_variable(model, abundance)
    community_biomass = com_model.metabolites.community_biomass_constraint
    abundance.add_metabolites({community_biomass: 1})

    print(i, com_model.optimize().f)
    com_model.repair()

    return com_model


def initialize_com_model(com_model, master_model, num_models, exchange_list):
    # Add exchanges from shared into individual model compartments
    # These need to be irreversible

    # Add exchanges for metabolites into shared community compartment
    for rxn in exchange_list:
        new_rxn = cobra.Reaction(rxn + '_shared')
        com_model.add_reaction(new_rxn)

        r = master_model.reactions.get_by_id(rxn)
        new_met = cobra.Metabolite(
            rxn.replace('EX_', '').replace('_reverse', '') + '_shared')
        for met, coeff in r.metabolites.items():
            new_rxn.add_metabolites({new_met: coeff})

        new_rxn.lower_bound = -1000
        new_rxn.upper_bound = 1000

    # Add variable for community biomass
    community_biomass_variable = cobra.Reaction('community_biomass_variable')
    community_biomass_variable.upper_bound = 1.
    com_model.add_reaction(community_biomass_variable)
    community_biomass_constraint = \
        cobra.Metabolite('community_biomass_constraint')
    community_biomass_variable.add_metabolites({community_biomass_constraint:
                                                    -1})


def run():
    rename_dict = {'lysA': 'b2838', 'metA': 'b4013', 'yddG': 'b1473',
                   'argH': 'b3960', 'pheA': 'b2599', 'yjeH': 'b4141',
                   'lysO': 'b0874', 'argO': 'b2923'}
    kos = [['lysA', 'metA', 'yddG'], ['argH', 'pheA', 'yjeH'],
           ['argH', 'lysO', 'pheA'], ['argO', 'lysA', 'metA']]

    kos = [['HISTD'], ['IPMD']]

    ijo = cobra.io.load_json_model(
        '/home/sbrg-cjlloyd/Desktop/ecoli_M_models/iJO1366.json')

    # add methionine export
    r = cobra.Reaction('METtpp')
    ijo.add_reaction(r)
    r.add_metabolites({'met__L_c': -1, 'met__L_p': 1})

    ijo.reactions.LYSt3pp.gene_reaction_rule = 'b0874'
    ijo.reactions.METtpp.gene_reaction_rule = 'b4141'

    com_model = cobra.Model('community_model')
    exchange_list = [i.id for i in ijo.reactions.query('EX_')]
    initialize_com_model(com_model, ijo, len(kos), exchange_list)

    model_list = [make_auxotrophs(ijo, ko, rename_dict) for ko in kos]
    for i, model in enumerate(model_list):
        com_model = add_model_to_community(com_model, model, i)
    return com_model
