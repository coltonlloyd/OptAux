from cobra import Model, Reaction, Metabolite
from cobra.manipulation.modify import convert_to_irreversible

from six import iteritems
from copy import deepcopy


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


def correct_equalities(model):
    # convert G and E constraints to L constraints
    for metabolite in model.metabolites:
        if metabolite._constraint_sense == "G":
            metabolite._constraint_sense = "L"
            metabolite._bound = - metabolite._bound
            for reaction in metabolite.reactions:
                coeff = reaction.get_coefficient(metabolite)
                # reverse the coefficient
                reaction.add_metabolites({metabolite: -2 * coeff})
        elif metabolite._constraint_sense == "E":
            # change existing constraint to L
            metabolite._constraint_sense = "L"
            # add new constraint
            new_constr = Metabolite("%s__GE_constraint" % metabolite.id)
            new_constr._constraint_sense = "L"
            new_constr._bound = - metabolite._bound
            for reaction in metabolite.reactions:
                coeff = reaction.get_coefficient(metabolite)
                reaction.add_metabolites({new_constr: -coeff})


def dual_problem(model, objective_sense="maximize", already_canonical=False,
                 integer_vars_to_maintain=[], even_dual=False,
                 already_irreversible=False, copy=True, dual_maximum=1000,
                 add_auxiliary=True):
    """Return a new model representing the dual of the model.

    Make the problem irreversible, then take the dual. Convert the problem:

    .. code-block:: none

        Maximize (c^T)x subject to Ax <= b, x >= 0

    which is something like this in COBRApy:

    .. code-block:: none

        Maximize sum(objective_coefficient_j * reaction_j for all j)
            s.t.
            sum(coefficient_i_j * reaction_j for all j) <= metabolite_bound_i
            reaction_j <= upper_bound_j
            reaction_j >= 0

    to the problem:

    .. code-block:: none

        Minimize (b^T)w subject to (A^T)w >= c, w >= 0

    which is something like this in COBRApy (S matrix is m x n):

    .. code-block:: none

        Minimize sum( metabolite_bound_i * dual_i   for all i ) +
                 sum( upper_bound_j *      dual_m+j for all j ) +
            s.t.
             sum( coefficient_i_j * dual_i for all i ) +
             sum( dual_2m+j' for all j' ) >= objective_coefficient_j
            dual_k >= 0


    Arguments
    ---------

    model : :class:`~cobra.core.Model`
        The COBRA model.

    objective_sense: str
        The objective sense of the starting problem, either 'maximize' or
        'minimize'. A minimization problems will be converted to a maximization
        before taking the dual. This function always returns a minimization
        problem.

    iteger_vars_to_maintain: [str]
        A list of IDs for Boolean integer variables to be maintained in the
        dual problem. See 'Maintaining integer variables' below for more
        details.

    already_irreversible: bool
        If True, then do not convert the model to irreversible.

    copy: bool
        If True, then make a copy of the model before modifying it. This is not
        necessary if already_irreversible is True.

    dual_maximum: float or int
        The upper bound for dual variables.


    **Maintaining integer variables**

    The argument ``integer_vars_to_maintain`` can be used to specify certin
    Boolean integer variables that will be maintained in the dual problem. This
    makes it possible to join outer and inner problems in a bi-level MILP. The
    method for maintaining integer variables is described by Tepper and Shlomi,
    2010:

    Tepper N, Shlomi T. Predicting metabolic engineering knockout strategies
    for chemical production: accounting for competing pathways. Bioinformatics.
    2010;26(4):536-43. https://doi.org/10.1093/bioinformatics/btp704.

    In COBRApy, this roughly translates to transforming (decision variables p,
    integer constraints o):

    .. code-block:: none

        Maximize (c^T)x subject to (A_x)x + (A_y)y <= b, x >= 0

        (1) Maximize sum(objective_coefficient_j * reaction_j for all j)
                s.t.
        (2)     sum(coeff_i_j * reaction_j for all j) +
                sum(decision_coeff_i_j * decision_var_j for all j)
                <= metabolite_bound_i
        (3)     reaction_j <= upper_bound_j
        (4)     reaction_j >= 0

    to the problem:

    .. code-block:: none

        Minimize (b - (A_y)y)^T w subject to (A_x^T)w >= c, w >= 0

    which linearizes to (with auxiliary variables z):

    .. code-block:: none

        Minimize (b^T)w - { ((A_y)y)^T w with yw --> z }
        subject to (A_x^T)w >= c, linearization constraints, w >= 0
          Linearization constraints: z <= w_max * y, z <= w,
                                     z >= w - w_max * (1 - y), z >= 0

        (5) Minimize sum( metabolite_bound_i *  dual_i            for all i ) +
                      sum( upper_bound_j *      dual_m+j          for all j ) +
                    - sum( decision_coeff_i_j * auxiliary_var_i_j
                          for all combinations i, j )
                s.t.
        (6)   - sum( coefficient_i_j * dual_i for all i ) - dual_m+j
              <= - objective_coefficient_j
        (7)     auxiliary_var_i_j - dual_maximum * decision_var_j          <= 0
        (8)     auxiliary_var_i_j - dual_i                                 <= 0
        (9)   - auxiliary_var_i_j + dual_i + dual_maximum * decision_var_j
              <= dual_maximum
       (10)     dual_maximum >= dual_i            >= 0
       (11)     dual_maximum >= dual_m+j          >= 0
       (12)     dual_maximum >= auxiliary_var_i_j >= 0
       (13)                1 >= decision_var_j    >= 0


    Zachary King 2015

    """

    # On every even dual transformation, we are taking the dual of a
    # minimization problem this is accounted for by changing the sign
    sign = -1 if even_dual else 1

    correct_equalities(model)

    # new model for the dual
    dual = Model("%s_dual" % model.id)

    # keep track of dual_i
    dual_var_for_met = {}

    # add dual variables for constraints. (2) --> dual_i
    for metabolite in model.metabolites:
        # add constraints based on metabolite constraint sense
        if metabolite._constraint_sense != "L":
            raise Exception("Not a less than or equal constraint: %s"
                            % metabolite.id)

        var = Reaction("%s__dual" % metabolite.id)
        # Without auxiliary variables, the objective coefficient would include
        # integer variables when present. However, we will separate out the
        # integer parts into objective coefficients for auxiliary variables.
        var.objective_coefficient = sign * metabolite._bound  # (5)
        # [dual_vars] >= 0
        # Set bounds
        var.lower_bound = 0
        var.upper_bound = dual_maximum

        dual.add_reaction(var)
        # remember
        dual_var_for_met[metabolite.id] = var

    # keep track of decision variables (integer_vars_to_maintain) as tuples:
    # (reaction in dual problem, reaction in original problem)
    integer_vars_added = []

    # add constraints and upper bound variables
    for reaction in model.reactions:

        # integer vars to maintain
        if reaction.id in integer_vars_to_maintain:
            # keep these integer variables in the dual, with new transformed
            # constraints
            if (reaction.lower_bound not in [0, 1] or
                    reaction.upper_bound not in [0, 1] or
                    reaction.variable_kind != "integer"):
                raise Exception("Reaction %s from integer_vars_to_maintain is "
                                "not a Boolean integer variable" % reaction.id)
            integer_var = Reaction(reaction.id)
            integer_var.upper_bound = reaction.upper_bound
            integer_var.lower_bound = reaction.lower_bound
            integer_var.variable_kind = reaction.variable_kind
            integer_var.objective_coefficient = 0
            # constraints
            dual.add_reaction(integer_var)
            integer_vars_added.append((integer_var, reaction))

        # other vars
        else:
            # other variables become constraints, (1) to (6)
            constr = Metabolite("%s__dual_constrained_by_c" %
                                reaction.id)  # (6)
            constr._constraint_sense = "G"
            # TODO make this negative to correct dual of dual problem
            constr._bound = sign * reaction.objective_coefficient
            for met, coeff in iteritems(reaction._metabolites):
                dual_var = dual_var_for_met[met.id]
                dual_var.add_metabolites({constr: coeff})

    # add auxiliary variables
    for integer_var, original_reaction in integer_vars_added:
        for metabolite, coeff in iteritems(original_reaction._metabolites):
            dual_var = dual_var_for_met[metabolite.id]
            # create an auxiliary variable
            aux_var = Reaction("%s__auxiliary__%s" % (integer_var.id,
                                                  dual_var.id))
            aux_var.lower_bound = 0
            aux_var.upper_bound = dual_maximum
            aux_var.objective_coefficient = - coeff
            dual.add_reaction(aux_var)

            # add linearization constraints
            # (7)     auxiliary_var_i_j - dual_maximum * decision_var_j    <= 0
            le_decision_constr = Metabolite("%s__le_decision" % aux_var.id)
            le_decision_constr._constraint_sense = "L"
            le_decision_constr._bound = 0
            aux_var.add_metabolites({le_decision_constr: 1})
            integer_var.add_metabolites({le_decision_constr: - dual_maximum})

            # (8)     auxiliary_var_i_j - dual_i                           <= 0
            le_dual_constr = Metabolite("%s__le_dual" % aux_var.id)
            le_dual_constr._constraint_sense = "L"
            le_dual_constr._bound = 0
            aux_var.add_metabolites({le_dual_constr: 1})
            dual_var.add_metabolites({le_dual_constr: -1})

            # (9)   - auxiliary_var_i_j + dual_i +
            #         dual_maximum * decision_var_j <= dual_maximum
            g_constr = Metabolite("%s__g_dual" % aux_var.id)
            g_constr._constraint_sense = "L"
            g_constr._bound = dual_maximum
            aux_var.add_metabolites({g_constr: -1})
            dual_var.add_metabolites({g_constr: 1})
            integer_var.add_metabolites({g_constr: dual_maximum})

    correct_equalities(dual)

    return dual


def canonical_form(model, objective_sense='maximize',
                   already_irreversible=False, copy=True):
    """Return a model (problem in canonical_form).

    Converts a minimization problem to a maximization, makes all variables
    positive by making reactions irreversible, and converts all constraints to
    <= constraints.


    model: class:`~cobra.core.Model`. The model/problem to convert.

    objective_sense: str. The objective sense of the starting problem, either
    'maximize' or 'minimize'. A minimization problems will be converted to a
    maximization.

    already_irreversible: bool. If the model is already irreversible, then pass
    True.

    copy: bool. Copy the model before making any modifications.

    """
    if copy:
        model = model.copy()

    if not already_irreversible:
        convert_to_irreversible(model)

    if objective_sense == "minimize":
        # if converting min to max, reverse all the objective coefficients
        for reaction in model.reactions:
            reaction.objective_coefficient = - reaction.objective_coefficient
    elif objective_sense != "maximize":
        raise Exception("Invalid objective sense '%s'. "
                        "Must be 'minimize' or 'maximize'." % objective_sense)

    correct_equalities(model)

    # convert lower bounds to LE constraints
    for reaction in model.reactions:
        # skip decision variables TODO fix this
        if reaction.variable_kind == 'integer':
            continue
        if reaction.lower_bound < 0:
            raise Exception("Bounds of irreversible reactions should be >= 0,"
                            " for %s" % reaction.id)
        # new constraint for lower bound
        lb_constr = Metabolite("%s__LB_constraint" % reaction.id)
        lb_constr._constraint_sense = "L"
        lb_constr._bound = - reaction.lower_bound
        reaction.add_metabolites({lb_constr: -1})

        # new constraint for lower bound
        lb_constr = Metabolite("%s__UB_constraint" % reaction.id)
        lb_constr._constraint_sense = "L"
        lb_constr._bound = reaction.upper_bound
        reaction.add_metabolites({lb_constr: 1})
        reaction.lower_bound = 0

    return model