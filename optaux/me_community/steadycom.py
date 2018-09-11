import cobra
import numpy as np
from cobrame import mu
from six import string_types, iteritems
from warnings import warn
from math import log
from cobra.solvers import solver_dict
from optaux.resources.possible_uptake import ko_uptakes
from cobra.flux_analysis.variability import calculate_lp_variability
from cobrame.solve.symbolic import substitute_mu, compile_expressions

from time import time

# Implementation of SteadyCom from Chan et al 2017
# (https://doi.org/10.1371/journal.pcbi.1005539)

try:
    import soplex
except ImportError as e:
    soplex = None
    warn("soplex import failed with error '%s'" % e.message)

try:
    from IPython.utils.coloransi import TermColors
    Red = TermColors.Red
    Green = TermColors.Green
    Normal = TermColors.Normal
except ImportError:
    Red = Green = Normal = ""


def make_auxotrophs(model_cons, ko, rename_dict= {}, ko_uptakes=[], ko_secretion=[]):
    model = model_cons.copy()
    if len(ko) > 0 and ko[0] in model_cons.reactions:
        for rxn in ko:
            model.reactions.get_by_id(rxn).knock_out()
    else:
        ko_bnums = [rename_dict[i] for i in ko]
        cobra.manipulation.remove_genes(model, ko_bnums)

    uptake_rxns = ko_uptakes
    if len(uptake_rxns) > 0:
        for r in uptake_rxns:
            model.reactions.get_by_id(r).lower_bound = -10.
    else:
        for r in model.reactions.query('EX_'):
            if r.lower_bound == 0:
                r.lower_bound = -10.

    if len(ko_secretion) > 0:
        for r in ko_secretion:
            model.reactions.get_by_id(r).upper_bound = 0

    model.id = str(ko)
    model.repair()
    print(model.id, model.optimize())
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

    for rxn in model.reactions:
        if rxn.id.startswith('EX_'):
            shared_met = rxn.id.replace('EX_', '') + '_shared'
            rxn.add_metabolites({cobra.Metabolite(shared_met): 1.})
            if rxn.lower_bound == 0:
                rxn.lower_bound = 0

    print(i, model.id, model.optimize().f)

    _add_explicit_bound_constraint(model)
    _add_constraints_to_abundance_strain_variable(model, abundance)

    # add metabolites to new model with appended ID for strain compartment
    for met in model.metabolites:
        if '_shared' not in met.id:
            met.id += '_strain_%i' % i
    model.repair()

    for rxn in model.reactions:
        new_rxn = rxn.copy()
        new_rxn.id += '_strain_%i' % i
        com_model.add_reaction(new_rxn)

    abundance = com_model.reactions.get_by_id('abundance_strain_%i' % i)

    community_biomass = com_model.metabolites.community_biomass_constraint
    abundance.add_metabolites({community_biomass: 1})

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

        new_rxn.lower_bound = -1000 if r.lower_bound < 0 else 0
        new_rxn.upper_bound = 1000

    # Add variable for community biomass
    community_biomass_variable = cobra.Reaction('community_biomass_variable')
    community_biomass_variable.upper_bound = 1.
    com_model.add_reaction(community_biomass_variable)
    community_biomass_constraint = \
        cobra.Metabolite('community_biomass_constraint')
    community_biomass_variable.add_metabolites({community_biomass_constraint:
                                                -1})


def get_ME_solver(solver=None):
    if solver is None:
        if soplex is None:
            raise RuntimeError("soplex not installed")
        return soplex
    elif isinstance(solver, string_types):
        return solver_dict[solver]
    else:
        return solver


def solve(model, min_mu=0, max_mu=2, mu_accuracy=1e-4,
          solver=None, verbose=True, compiled_expressions=None,
          debug=True, reset_obj=False, **solver_args):
    """Modified from COBRAme"""

    if solver is not None:
        debug = False  # other solvers can't handle debug mode
    solver = get_ME_solver(solver)
    lp = solver.create_problem(model)
    # reset the objective for faster feasibility solving
    if reset_obj:
        for i, reaction in enumerate(model.reactions):
            if reaction.objective_coefficient != 0:
                solver.change_variable_objective(lp, i, 0)
    for name, value in iteritems(solver_args):
        solver.set_parameter(lp, name, value)
    if compiled_expressions is None:
        compiled_expressions = compile_expressions(model)
    feasible_mu = []
    infeasible_mu = []

    # String formatting for display
    str_places = int(abs(round(log(mu_accuracy) / log(10)))) + 1
    num_format = "%." + str(str_places) + "f"
    mu_str = "mu".ljust(str_places + 2)
    if debug:
        verbose = True
    if verbose:
        success_str_base = Green + num_format + "\t+" + Normal
        failure_str_base = Red + num_format + "\t-" + Normal
        if debug:
            print("%s\tstatus\treset\ttime\titer\tobj" % mu_str)
        else:
            print("%s\tstatus" % mu_str)

    def try_mu(mu):
        substitute_mu(lp, mu, compiled_expressions, solver)

        solver.solve_problem(lp)
        status = solver.get_status(lp)
        print(status)
        if debug:
            reset_basis = lp.reset_basis
            obj = str(lp.get_objective_value()) \
                if status == "optimal" else ""
            debug_str = "\t%s\t%.2f\t%d\t%s" % \
                        (reset_basis, lp.solveTime, lp.numIterations, obj)
        else:
            debug_str = ""
        if status == 'optimal':
            if verbose:
                print(success_str_base % mu + debug_str)
            feasible_mu.append(mu)
            return True
        else:
            print(failure_str_base % mu + debug_str)
            infeasible_mu.append(mu)
            return False

    start = time()
    # find highest possible value of mu try the edges of binary search
    if not try_mu(min_mu):
        # Try 0 if min_mu failed
        if min_mu == 0 or not try_mu(0):
            raise ValueError("0 needs to be feasible")
    while try_mu(max_mu):  # If max_mu was feasible, keep increasing
        max_mu += 1
    while infeasible_mu[-1] - feasible_mu[-1] > mu_accuracy:
        try_mu((infeasible_mu[-1] + feasible_mu[-1]) * 0.5)
    # now we want to solve with the objective
    if reset_obj:
        for i, reaction in enumerate(model.reactions):
            if reaction.objective_coefficient != 0:
                solver.change_variable_objective(
                    lp, i, reaction.objective_coefficient)
    try_mu(feasible_mu[-1])
    model.solution = solver.format_solution(lp, model)
    model.solution.f = feasible_mu[-1]

    if verbose:
        print("completed in %.1f seconds and %d iterations" %
              (time() - start, len(feasible_mu) + len(infeasible_mu)))


def run_fva_for_abundances(model, percent_max, max_gr, solver=None,
                           compiled_expressions=None):
    for strain in model.reactions.query('abundance_strain_'):
        strain.lower_bound = 0
        strain.upper_bound = 1.

    for r in model.reactions:
        r.objective_coefficient = 0

    solver = get_ME_solver(solver)
    lp = solver.create_problem(model)
    # reset the objective for faster feasibility solving

    if compiled_expressions is None:
        compiled_expressions = compile_expressions(model)
    from collections import defaultdict
    out_abundances = defaultdict(dict)
    for strain in model.reactions.query('abundance_strain_'):
        for frac in np.linspace(percent_max, 1, 10):
            substitute_mu(lp, frac*max_gr, compiled_expressions, solver)
            out_dict = calculate_lp_variability(lp, solver, model, [strain])
            out_abundances[frac].update(out_dict)

    return out_abundances


def run(ko1_list, ko2_list):

    kos = [ko1_list, ko2_list]

    ijo = cobra.io.load_json_model(
        '/home/sbrg-cjlloyd/Desktop/ecoli_M_models/iJO1366.json')

    ijo.reactions.EX_o2_e.lower_bound = -1000
    ijo.reactions.ATPM.lower_bound = 3.15
    ijo.reactions.EX_glc__D_e.lower_bound = -10
    # add methionine export
    r = cobra.Reaction('OROTtpp')
    ijo.add_reaction(r)
    r.add_metabolites({'orot_c': -1, 'orot_p': 1})

    com_model = cobra.Model('community_model')
    exchange_list = [i.id for i in ijo.reactions.query('EX_')]
    initialize_com_model(com_model, ijo, len(kos), exchange_list)

    model_list = [make_auxotrophs(ijo, ko, ko_uptakes=ko_uptakes[str(ko)]) for i, ko in enumerate(kos)]
    for i, model in enumerate(model_list):
        com_model = add_model_to_community(com_model, model, i)

    com_model.reactions.community_biomass_variable.lower_bound = 1.
    return com_model
    # solve(com_model)


def reproduce_figure():
    rename_dict = {'lysA': 'b2838', 'metA': 'b4013', 'yddG': 'b1473',
                   'argH': 'b3960', 'pheA': 'b2599', 'yjeH': 'b4141',
                   'lysO': 'b0874', 'argO': 'b2923'}

    ijo = cobra.io.load_json_model(
        '/home/sbrg-cjlloyd/Desktop/ecoli_M_models/iJO1366.json')
    # add methionine export
    r = cobra.Reaction('METtpp')
    ijo.add_reaction(r)
    r.add_metabolites({'met__L_c': -1, 'met__L_p': 1})
    ijo.reactions.LYSt3pp.gene_reaction_rule = 'b0874'
    ijo.reactions.METtpp.gene_reaction_rule = 'b4141'

    kos = [['lysA', 'metA', 'yddG'],
           ['argH', 'pheA', 'yjeH'],
           ['argH', 'lysO', 'pheA'],
           ['argO', 'lysA', 'metA']]

    possible_uptake = [['EX_lys__L_e', 'EX_met__L_e'],
                       ['EX_arg__L_e', 'EX_phe__L_e'],
                       ['EX_arg__L_e', 'EX_phe__L_e'],
                       ['EX_lys__L_e', 'EX_met__L_e']]

    # metabolites each strain cannot secrete
    ko_secretion = [['EX_lys__L_e', 'EX_met__L_e', 'EX_phe__L_e'],
                    ['EX_arg__L_e', 'EX_met__L_e', 'EX_phe__L_e'],
                    ['EX_arg__L_e', 'EX_lys__L_e', 'EX_phe__L_e'],
                    ['EX_arg__L_e', 'EX_lys__L_e', 'EX_met__L_e']]
    ijo.reactions.EX_o2_e.lower_bound = -18.5
    ijo.reactions.ATPM.lower_bound = 10.39
    ijo.reactions.EX_glc__D_e.lower_bound = -8

    com_model = cobra.Model('community_model')
    exchange_list = [i.id for i in ijo.reactions.query('EX_')]
    initialize_com_model(com_model, ijo, len(kos), exchange_list)

    model_list = [make_auxotrophs(ijo, ko, ko_uptakes=possible_uptake[i],
                                  ko_secretion=ko_secretion[i],
                                  rename_dict=rename_dict) for
                  i, ko in enumerate(kos)]
    for i, model in enumerate(model_list):
        com_model = add_model_to_community(com_model, model, i)

    com_model.reactions.community_biomass_variable.lower_bound = 1.
    return com_model