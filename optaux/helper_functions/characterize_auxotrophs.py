import cobra
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
import numpy as np


def get_auxotrophic_mets_per_ko(cons_model, KOs, growth_rate=.1):
    metabolite_list = []
    model = cons_model.copy()
    for r in KOs:
        model.reactions.get_by_id(r).knock_out()

    if model.optimize().f > growth_rate:
        return []

    solver = cobra.solvers.solver_dict['gurobi']
    lp = solver.create_problem(model)
    for rxn in model.reactions.query('EX_'):
        old_bounds = (rxn.lower_bound, rxn.upper_bound)
        index = model.reactions.index(rxn)
        solver.change_variable_bounds(lp, index, -10., old_bounds[1])
        solver.solve_problem(lp)
        # get the status and growth rate
        status = solver.get_status(lp)
        # reset the problem
        solver.change_variable_bounds(lp, index, old_bounds[0], old_bounds[1])
        f = solver.get_objective_value(lp) if status == "optimal" else 0.
        if f > .1:
            metabolite_list.append(rxn.id)
    return metabolite_list


def get_avg_flux_required(cons_model, KOs, aux_metabolite_list):
    fluxes = []

    model = cons_model.copy()
    biomass = list(model.objective.keys())[0]
    for r2 in KOs:
        model.reactions.get_by_id(r2).knock_out()
    for r in aux_metabolite_list:
        if r == 'EX_h2o_e':
            continue
        r_obj = model.reactions.get_by_id(r)
        model.objective = biomass
        r_obj.lower_bound = -10

#        sol = model.optimize()
        biomass.lower_bound = .1
        biomass.upper_bound = .1
        model.objective = r_obj

        sol2 = model.optimize()
        try:
            print(r, sol2.x_dict[r])
            fluxes.append(sol2.x_dict[r])
        except:
            print(r, KOs, 'ISSUE')
        r_obj.lower_bound = 0.
        biomass.lower_bound = 0.
        biomass.upper_bound = 1000.

    if len(fluxes) > 0:
        return np.array(fluxes).mean()
    else:
        return None


def get_blocked_biomass(cons_model, KOs):
    blocked = []
    model = cons_model.copy()

    biomass = list(model.objective.keys())[0]
    for r in KOs:
        model.reactions.get_by_id(r).knock_out()
    for metabolite in biomass.reactants:
        demand = cobra.Reaction('DM_' + metabolite.id)
        demand.add_metabolites({metabolite: -1})
        model.add_reaction(demand)
        model.change_objective(demand)
        sol2 = model.optimize()
        if sol2.f < .1:
            blocked.append(metabolite.id)
    return blocked


def gene_names_per_kos(model, KOs):
    gene_names = []
    for ko in KOs:
        print(model.reactions.get_by_id(ko).gene_name_reaction_rule)
        gene_names.append(
            model.reactions.get_by_id(ko).gene_name_reaction_rule)
    return gene_names


def genes_per_kos(model, KOs):
    genes = []
    for ko in KOs:
        genes.append(model.reactions.get_by_id(ko).gene_reaction_rule)
    return genes


def to_string(list):
    return ', '.join(list)


def and_join_strings(list):
    return ' & '.join(list)


def id_to_name(model, ids):
    names = []
    for i in ids:
        if i.startswith('EX_'):
            names.append(model.reactions.get_by_id(i).name)
        else:
            names.append(model.metabolites.get_by_id(i).name)
    return names
