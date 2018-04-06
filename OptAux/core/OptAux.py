from __future__ import print_function, absolute_import

from copy import deepcopy
from OptAux.core.dual import dual_problem, canonical_form
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
        lower_bound = - model.reactions.get_by_id(reaction_id + '_reverse').upper_bound
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


def set_up_optaux_old(model, chemical_objective, knockable_reactions,
                  min_biomass=.1, n_knockouts=1, dual_maximum=1000.,
                  n_knockouts_required=False, type='max-min'):
    raise DeprecationWarning("don't use this")
    model = model.copy()

    decision_variable_ids = [
        design.design_algorithms._add_decision_variable(model, r_id).id
        for r_id in knockable_reactions]

    objective_rxn = model.reactions.get_by_id(chemical_objective)

    model.reactions.BIOMASS_Ec_iJO1366_core_53p95M.lower_bound = min_biomass
    model.reactions.BIOMASS_Ec_iJO1366_core_53p95M.upper_bound = min_biomass

    print(model.objective)

    # inner problem
    inner_problem = model.copy()

    if objective_rxn.objective_coefficient == 1:
        objective_sense = 'maximize'
    elif objective_rxn.objective_coefficient == -1:
        objective_sense = 'minimize'
    else:
        raise('Must have chemical objective')

    inner_dual = design.dual_problem(inner_problem,
                                     integer_vars_to_maintain=decision_variable_ids,
                                     already_irreversible=False, copy=False,
                                     dual_maximum=dual_maximum,
                                     objective_sense='minimize')

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

    # constraint to set outer and inner objectives equal, and set chemical
    # objective
    equal_objectives_constr = Metabolite("equal_objectives_constraint")
    equal_objectives_constr._constraint_sense = "E"
    equal_objectives_constr._bound = 0
    for reaction in model.reactions:
        if reaction.objective_coefficient != 0:
            reaction.add_metabolites({equal_objectives_constr:
                                          - reaction.objective_coefficient})
        inner_objective = inner_objectives.get(reaction.id, 0)
        if inner_objective:
            reaction.add_metabolites(
                {equal_objectives_constr: - inner_objective})
        # set chemical objective
        #reaction.objective_coefficient = -1 \
        #    if reaction.id == chemical_objective else 0

    # add the n_knockouts constraint
    n_knockouts_constr = Metabolite("n_knockouts_constraint")
    n_knockouts_constr._constraint_sense = "E" if n_knockouts_required else "G"
    n_knockouts_constr._bound = len(decision_variable_ids) - n_knockouts
    for r_id in decision_variable_ids:
        reaction = model.reactions.get_by_id(r_id)
        reaction.add_metabolites({n_knockouts_constr: 1})

    return model


def dual_embed(cons_model, decision_variable_ids, embed_type='max-max',
               even_dual=False, dual_maximum=1000., already_canonical=True):
    """Possible inner_types = 'max-max' or 'max-min'"""

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

    ## constraint to set outer and inner objectives equal
    #equal_objectives_constr = Metabolite("equal_objectives_constraint_1")
    #equal_objectives_constr._constraint_sense = "L"
    #equal_objectives_constr._bound = 0
    #for reaction in model.reactions:
    #    if reaction.objective_coefficient != 0:
#
    #        reaction.add_metabolites({equal_objectives_constr:
    #                                      - reaction.objective_coefficient})
#
    #    inner_objective = inner_objectives.get(reaction.id, 0)
    #    if inner_objective:
    #        reaction.add_metabolites(
    #            {equal_objectives_constr: inner_objective})
#
    return model


def set_up_design_problem(cons_model, chemical_objective, knockable_reactions,
                          min_biomass=.1, n_knockouts=1, dual_maximum=1000.,
                          n_knockouts_required=True,
                          optimization_type='OptKnock'):
    model = cons_model.copy()

    decision_variable_ids = [_add_decision_variable(model, r_id).id
                             for r_id in knockable_reactions]

    for r in model.reactions:
        if r.objective_coefficient != 0:
            growth_objective = r.id
            r.objective_coefficient = 0

    model.reactions.get_by_id(growth_objective).lower_bound = min_biomass
    model.reactions.get_by_id(growth_objective).upper_bound = 1000

    print(model.objective)

    # convert to canonical form and copy
    con_model = canonical_form(model, copy=True)

    cons_model = con_model.copy()

    con_model.change_objective(growth_objective)

    model = dual_embed(con_model, decision_variable_ids,
                       embed_type='max-max', even_dual=False,
                       dual_maximum=dual_maximum, already_canonical=True)

    for r in model.reactions:
        if r.id == chemical_objective:
            r.objective_coefficient = 1
        elif r.id == chemical_objective + '_reverse':
            r.objective_coefficient = -1
        else:
            r.objective_coefficient = 0.

    if optimization_type == 'RobustKnock':

        for r in model.reactions:
            if r.id == chemical_objective:
                r.objective_coefficient = -1
            elif r.id == chemical_objective + '_reverse':
                r.objective_coefficient = 1
            else:
                r.objective_coefficient = 0.

        model = dual_problem(model,
                             integer_vars_to_maintain=decision_variable_ids)

        for r in cons_model.reactions:
            if 'decision_var' not in r.id:
                r.objective_coefficient = 0
                model.add_reaction(r)
            else:
                r.objective_coefficient = 0
        for r in cons_model.reactions:
            if 'decision_var' in r.id:
                update_decision_variable(model,
                                         r.id.replace('_decision_var', ''))

        # max c_dual^{-T} * x_dual
        for r in model.reactions:
            r.objective_coefficient = -r.objective_coefficient

    # add the n_knockouts constraint
    n_knockouts_constr = Metabolite("n_knockouts_constraint")
    n_knockouts_constr._constraint_sense = "E" if n_knockouts_required else "G"
    n_knockouts_constr._bound = len(decision_variable_ids) - n_knockouts
    for r_id in decision_variable_ids:
        reaction = model.reactions.get_by_id(r_id)
        reaction.add_metabolites({n_knockouts_constr: 1})
    return model, decision_variable_ids


def set_up_optaux(model, chemical_objective, knockable_reactions,
                  min_biomass=.1, n_knockouts=1, dual_maximum=1000.,
                  n_knockouts_required=False, type='max-min'):
    """Use this implementation. It is more correct."""
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

    for r in con_model.reactions:
        if r.id == chemical_objective:
            r.objective_coefficient = 1
        elif r.id == chemical_objective + '_reverse':
            r.objective_coefficient = -1
        else:
            r.objective_coefficient = 0.

    model = con_model.copy()
    cons_model = con_model.copy()

    model = dual_problem(model, integer_vars_to_maintain=decision_variable_ids)
    for r in cons_model.reactions:
        if 'decision_var' not in r.id:
            r.objective_coefficient = 0
            model.add_reaction(r)
        else:
            r.objective_coefficient = 0
    for r in cons_model.reactions:
        if 'decision_var' in r.id:
            update_decision_variable(model, r.id.replace('_decision_var', ''))

    # max c_dual^{-T} * x_dual
    for r in model.reactions:
        r.objective_coefficient = -r.objective_coefficient

    # add the n_knockouts constraint
    n_knockouts_constr = Metabolite("n_knockouts_constraint")
    n_knockouts_constr._constraint_sense = "E" if n_knockouts_required else "G"
    n_knockouts_constr._bound = len(decision_variable_ids) - n_knockouts
    for r_id in decision_variable_ids:
        reaction = model.reactions.get_by_id(r_id)
        reaction.add_metabolites({n_knockouts_constr: 1})

    return model
