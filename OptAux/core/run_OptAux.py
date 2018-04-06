from __future__ import print_function, division, absolute_import
import time

from OptAux.core.OptAux import set_up_optaux
from OptAux.core.ReduceModel import reduce_model, get_targets

subsystem_exclude_list = ['Cell Envelope Biosynthesis',
                          'Exchange',
                          'Inorganic Ion Transport and Metabolism',
                          'Lipopolysaccharide Biosynthesis / Recycling',
                          'Murein Biosynthesis',
                          'Murein Recycling',
                          'Transport, Inner Membrane',
                          'Transport, Outer Membrane',
                          'Transport, Outer Membrane Porin',
                          'tRNA Charging'
                          ]


def get_knockouts(optknock_problem):
    knockouts = []
    solution = optknock_problem.solution
    for reaction in optknock_problem.reactions:
        if solution.x_dict.get(reaction.id, None) == 0:
            r_id = getattr(reaction, "decision_reaction_id", None)
            if r_id is not None:
                knockouts.append(r_id)
    return knockouts


def run_optaux(cons_model, target_metabolite, percent_max_growth,
               n_knockouts, aerobicity='aerobic',
               lb_essential_list=[], trace_metabolite_threshold=0,
               media_list=[], exclude_reactions={},
               solver='gurobi'):

    model = cons_model.copy()

    # Find growth objective
    growth_objective = list(model.objective.keys())[0].id

    # Target exchange reaction
    target_exchange = 'EX_' + target_metabolite

    # Set oxygen uptake if aerobic
    if aerobicity == 'aerobic':
        model.reactions.EX_o2_e.lower_bound = -20.
        oxygen_bound = -20
        print('-20 oxygen')
    elif aerobicity == 'anaerobic':
        model.reactions.EX_o2_e.lower_bound = 0.
        oxygen_bound = 0
    else:
        raise ValueError('Must set aerobicity')

    if trace_metabolite_threshold:
        for r in model.reactions.query('EX_'):
            if r.id == 'EX_glc__D_e' or r.lower_bound < 0:
                continue
            r.lower_bound = -trace_metabolite_threshold

    # Scale target metabolite uptake by number of carbons
    scale = model.metabolites.glc__D_c.elements['C'] / \
            model.metabolites.get_by_id(target_metabolite).elements['C']

    # Add target metabolite to media
    model.reactions.get_by_id(target_exchange).lower_bound = \
        -10 #* scale

    # Add metabolites present in media, if any
    for met in media_list:
        scale = model.metabolites.glc__D_c.elements['C'] / \
                model.metabolites.get_by_id(met).elements['C']
        model.reactions.get_by_id('EX_' + met).lower_bound = -10 * scale

    sol = model.optimize(solver=solver)
    if sol.f < .5:
        raise Exception('Error preventing model from solving')

    # Reduce model and get set of possible reactions for knock out
    model_red = reduce_model(model, change_bounds=False, tol=1e-9)
    knockable_reactions = \
        get_targets(model_red, subsys_exclude_list=subsystem_exclude_list,
                    manual_exclude_reaction_list=lb_essential_list)

    # Set the min_biomass to a percent of the max growth rate
    max_WT_growth = model_red.optimize().f
    min_biomass = max_WT_growth * percent_max_growth

    # Set objective to uptake of target metabolite
    try:
        model_red.reactions.get_by_id(target_exchange).objective_coefficient = 0
    except:
        print('Reaction (%s) not in reduced model' % target_exchange)
        return {}
    model_red.reactions.get_by_id(growth_objective).objective_coefficient = 1.

    print(min_biomass)
    # Set up the final OptAux problem
    p = set_up_optaux(model_red, target_exchange,
                      set(knockable_reactions)-set(exclude_reactions),
                      n_knockouts=n_knockouts, min_biomass=min_biomass,
                      n_knockouts_required=True)

    # Optimize
    tic = time.time()
    print(p.optimize(solver=solver))
    toc = time.time() - tic

    knockouts = get_knockouts(p)
    # Test the OptAux solution to find the minimum required uptake
    model_test = model_red.copy()
    for r in knockouts:
        model_test.reactions.get_by_id(r).knock_out()

    growth_reaction = model_test.reactions.get_by_id(growth_objective)
    growth_reaction.lower_bound = min_biomass
    growth_reaction.upper_bound = min_biomass
    model_test.objective = target_exchange
    min_uptake = model_test.optimize('maximize').f

    # find max growth rate
    growth_reaction.upper_bound = 1000.
    growth_reaction.lower_bound = 0.
    model_test.objective = growth_objective
    gr = model_test.optimize().f

    # find min uptake at max growth rate
    model_test.objective = target_exchange
    growth_reaction.lower_bound = gr
    growth_reaction.upper_bound = gr
    min_max_uptake = model_test.optimize().f

    return {'Reaction Knockouts': knockouts,
            'Minimum Uptake at set_biomass':
                min_uptake if min_uptake < -1e-12 else 0,  # account for floating point precision
            'Number of Knockouts': n_knockouts,
            'Growth Rate used for set_biomass': min_biomass,
            'Maximum Growth Rate w/ Knockouts': gr,
            'Minimum Uptake at Max Growth Rate w/ Knockouts':
                min_max_uptake if min_max_uptake < -1e-12 else 0,
            'Oxygen Uptake': oxygen_bound,
            'Time (s)': toc,
            'Trace Metabolite Threshold': trace_metabolite_threshold,
            'Solver': solver, 'Solve Status': p.solution.status}
