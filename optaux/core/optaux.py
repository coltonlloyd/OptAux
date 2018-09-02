from __future__ import print_function, absolute_import

from copy import deepcopy
from optaux.core.dual import dual_problem, canonical_form
from cobra.core import Model, Reaction, Metabolite


def _add_decision_variable(model, reaction_id):
    """Add an integer decision variable for the given reaction."""
    reaction = model.reactions.get_by_id(reaction_id)
    # add integer variable
    var = Reaction("%s_decision_var" % reaction_id)
    var.lower_bound = 0
    var.upper_bound = 1
    var.variable_kind = "integer"
    var.decision_reaction_id = reaction_id
    model.add_reaction(var)
    # add constraints
    # v <= ub * y  -->  v - ub * y <= 0
    ub_constr = Metabolite("%s__UB_constraint" % var.id)
    ub_constr._constraint_sense = "L"
    # v >= lb * y  -->  - v + lb * y <= 0
    lb_constr = Metabolite("%s__LB_constraint" % var.id)
    lb_constr._constraint_sense = "L"
    reaction.add_metabolites({lb_constr: - 1,
                              ub_constr:   1})
    var.add_metabolites({lb_constr:   reaction.lower_bound,
                         ub_constr: - reaction.upper_bound})
    return var


def update_decision_variable(model, reaction_id):
    """Add an integer decision variable for the given reaction.

    Ensure decision variable applies to forward and reverse direction"""
    reaction = model.reactions.get_by_id(reaction_id)
    upper_bound = reaction.upper_bound

    if reaction_id + '_reverse' in model.reactions:
        lower_bound = \
            - model.reactions.get_by_id(reaction_id + '_reverse').upper_bound
    else:
        lower_bound = reaction.lower_bound
    # add integer variable
    var = model.reactions.get_by_id("%s_decision_var" % reaction_id)
    var.decision_reaction_id = reaction_id

    # add constraints
    # v <= ub * y  -->  v - ub * y <= 0
    ub_constr = Metabolite("%s__UB_constraint" % var.id)
    ub_constr._constraint_sense = "L"
    # v >= lb * y  -->  - v + lb * y <= 0
    lb_constr = Metabolite("%s__LB_constraint" % var.id)
    lb_constr._constraint_sense = "L"

    var.add_metabolites({lb_constr: lower_bound,
                         ub_constr: - upper_bound})
    return var


def dual_embed(cons_model, decision_variable_ids,
               even_dual=False, dual_maximum=1000., already_canonical=True):
    """Dual embed max-mix"""

    # inner problem
    inner_problem = cons_model.copy()
    model = cons_model.copy()

    inner_dual = dual_problem(inner_problem,
                              integer_vars_to_maintain=decision_variable_ids,
                              already_irreversible=False, copy=False,
                              dual_maximum=dual_maximum,
                              even_dual=even_dual,
                              already_canonical=already_canonical)

    # add constraints and variables from inner problem to outer problem
    inner_objectives = {}
    for reaction in inner_dual.reactions:
        inner_objectives[reaction.id] = reaction.objective_coefficient
        reaction.objective_coefficient = 0
        if reaction.id in model.reactions:
            existing_reaction = model.reactions.get_by_id(reaction.id)
            for met, coeff in reaction._metabolites.items():
                if met.id in model.metabolites:
                    existing_reaction.add_metabolites(
                        {model.metabolites.get_by_id(met.id): coeff})
                else:
                    existing_reaction.add_metabolites({deepcopy(met): coeff})
        else:
            model.add_reaction(reaction)

    # constraint to set outer and inner objectives equal
    equal_objectives_constr = Metabolite("equal_objectives_constraint")
    equal_objectives_constr._constraint_sense = "E"
    equal_objectives_constr._bound = 0
    for reaction in model.reactions:
        if reaction.objective_coefficient != 0:
            reaction.add_metabolites({equal_objectives_constr:
                                      reaction.objective_coefficient})
        inner_objective = inner_objectives.get(reaction.id, 0)
        if inner_objective:
            reaction.add_metabolites(
                {equal_objectives_constr: -inner_objective})

    return model


def set_up_optaux(model, chemical_objective, knockable_reactions,
                  min_biomass=.1, n_knockouts=1, dual_maximum=1000.,
                  n_knockouts_required=False):
    """Sets up OptAux problem for optimization."""
    model = model.copy()

    decision_variable_ids = [_add_decision_variable(model, r_id).id
                             for r_id in knockable_reactions]

    objective_rxn = model.reactions.get_by_id(chemical_objective)
    objective_rxn.lower_bound = -20
    print(model.objective)
    growth_rxn = list(model.objective.keys())[0]

    growth_rxn.lower_bound = min_biomass
    growth_rxn.upper_bound = min_biomass

    # convert to canonical form and copy
    con_model = canonical_form(model, copy=True)

    # Inner problem is a minimization
    for r in con_model.reactions:
        if r.id == chemical_objective:
            r.objective_coefficient = 1
        elif r.id == chemical_objective + '_reverse':
            r.objective_coefficient = -1
        else:
            r.objective_coefficient = 0.

    model = con_model.copy()

    # Dual embed minimization problem
    model = dual_embed(model, decision_variable_ids)

    # Reset out objective as maximization
    for r in model.reactions:
        if r.id == chemical_objective:
            r.objective_coefficient = -1
        elif r.id == chemical_objective + '_reverse':
            r.objective_coefficient = 1
        else:
            r.objective_coefficient = 0.

    # add the n_knockouts constraint
    n_knockouts_constr = Metabolite("n_knockouts_constraint")
    n_knockouts_constr._constraint_sense = "E" if n_knockouts_required else "G"
    n_knockouts_constr._bound = len(decision_variable_ids) - n_knockouts
    for r_id in decision_variable_ids:
        reaction = model.reactions.get_by_id(r_id)
        reaction.add_metabolites({n_knockouts_constr: 1})

    return model
