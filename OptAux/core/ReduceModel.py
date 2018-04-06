from __future__ import print_function
import cobra


currencyMets = {'h2o', 'co2', 'o2', 'h2o2', 'nh4', 'no2', 'no3', 'no', 'h2s',
                'so3','so4','h','h2','pi','ppi','coa','accoa','ppcoa','aacoa',
                'butcoa','succoa','atp','gtp','adp','gdp','amp','gmp','nad',
                'nadp','nadh','nadph','fad','fadh','na1','ahcys','amet','thf','mlthf',
                'q8h2','q8','mql8','mqn8','2dmmql8','2dmmq8'}


def check_consistency(model, model_red, tol, solver='gurobi'):
    diff_max = model_red.optimize(solver=solver).f - model.optimize(solver=solver).f
    diff_min = model_red.optimize('minimize', solver=solver).f - model.optimize('minimize', solver=solver).f
    consistency = True
    if abs(diff_max) > tol:
        consistency = False
        print('Max wrong %s, %s' % (model_red.optimize().f, model.optimize().f))
    if abs(diff_min) > tol:
        consistency = False
        print('Min wrong %s, %s' % (model_red.optimize('minimize').f,
                                    model.optimize('minimize').f))

    return consistency


def expand_bounds(model, model_red, tol):
    cushion = tol

    model_ok = False
    while not model_ok:
        for r in model_red.reactions:
            if (r.upper_bound - r.lower_bound) < cushion and \
                            r.upper_bound != r.lower_bound:
                if r.lower_bound != 0:
                    r.lower_bound = r.lower_bound - cushion
                if r.upper_bound != 0:
                    r.upper_bound = r.upper_bound + cushion
        cushion = cushion * 2.
        model_ok = check_consistency(model, model_red, tol)
    return model_red


def reduce_model(model_cons, tol=1e-8, irreversible_flag=False,
                 change_bounds=True, solver='gurobi'):
    objective = list(model_cons.objective.keys())[0].id
    model = model_cons.copy()
    model_red = model_cons.copy()
    # cobra.manipulation.convert_to_irreversible(model)
    # cobra.manipulation.convert_to_irreversible(model_red)
    cobra.io.save_json_model(model_red, 'model.json')
    fva_sol = cobra.flux_analysis.flux_variability_analysis(model, solver=solver,
                                                            fraction_of_optimum=0.)

    for rxn_id, solution in fva_sol.items():
        rxn = model_red.reactions.get_by_id(rxn_id)
        lower_bound = solution['minimum'] if abs(
            solution['minimum']) > tol else 0
        upper_bound = solution['maximum'] if abs(
            solution['maximum']) > tol else 0
        if upper_bound == 0. and lower_bound == 0.:
            rxn.remove_from_model(remove_orphans=False)
        elif change_bounds:
            rxn.lower_bound = lower_bound
            rxn.upper_bound = upper_bound

    # Check consistency
    model.objective = objective
    model_red.objective = objective

    cobra.io.save_json_model(model_red, 'model_red.json')

    modelOK = check_consistency(model, model_red, tol)

    if not modelOK:
        model_red = expand_bounds(model, model_red, tol)

    return model_red


def get_targets(model, subsys_exclude_list=[], num_carbons_threshold=10,
                manual_exclude_reaction_list=[]):
    max_gr = model.optimize().f

    spont = model.genes.get_by_id('s0001')
    possible_KO_reactions = []
    for r in model.reactions:
        if r.id.startswith('DM_') or r.id.startswith('EX_'):
            continue
        if r.subsystem in subsys_exclude_list:
            continue
        if r.id in manual_exclude_reaction_list:
            continue
        if spont in r.genes or r.gene_reaction_rule == '':
            continue

        high_carbon = False
        for met in r.metabolites:
            if met.id[:-2] in currencyMets:
                continue
            if met.elements.get('C', 0) > num_carbons_threshold:
                high_carbon = True
        if high_carbon:
            continue
        possible_KO_reactions.append(r)

    ko_gr = cobra.flux_analysis.single_reaction_deletion(model,
                                                         possible_KO_reactions)

    knockable_reactions = [i for i, v in ko_gr[0].items() if v > max_gr * .05]

    return knockable_reactions
