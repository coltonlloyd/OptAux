# ==========================================================
# Sensitivity of substrate uptake hierarchy to ME params
#
# Laurence Yang, SBRG, UCSD
#
# 25 Mar 2016: first version
# 08 Apr 2016: sample starting at calibrated model (keffs)
# 17 Aug 2016: modified by Colton Lloyd

# =============== Supported simulation types =================
# carbon_substrates: all growth supporting carbon sources

# ============================================================
from qminospy.me1 import ME_NLP1
import numpy as np
import pickle
import os
from os.path import abspath, dirname, relpath
import json
import time
import pandas as pd
import argparse

# ------------------------------------------------------------
# Modules to Manipulate Model
# ------------------------------------------------------------
from optaux.me_community.run_ALE_pairs import setup_simulation
from optaux import resources
from optaux.resources.possible_uptake import get_possible_uptake
import cobra
# ************************************************************
# Parameters
# ------------------------------------------------------------
# Bisection parameters
MU_PREC = 1e-4
MU_MIN = 0.
MU_MAX = .7

default_secretion_keff = 6.5

simstr = ''


def str2bool(input):
    if input.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif input.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# ************************************************************
# Argument
parser = argparse.ArgumentParser(description='Simulation parameters.')

parser.add_argument('pair', help='', type=str)
parser.add_argument('fraction', help='fraction of strain 1', type=float)
parser.add_argument('unmodeled1', help='unmodeled protein fraction strain 1',
                    type=float)
parser.add_argument('unmodeled2', help='unmodeled protein fraction strain 2',
                    type=float)
parser.add_argument('mode', help='the type of simulations being ran',
                    type=str)
parser.add_argument('--keff_transporter_1', help='multiplier of outer membrane'
                                                 'tansporter of strain 1',
                    type=float)
parser.add_argument('--keff_transporter_2', help='multiplier of outer membrane'
                                                 'tansporter of strain 2',
                    type=float)
parser.add_argument('--Scale_secretion', default='True')
parser.add_argument('--Restrict_crossfeeding', default='experimental_inferred')
parser.add_argument('--glucose_uptake', default=-1000)
parser.add_argument('--docker', default='False')

args = parser.parse_args()

# pair must be in form of ko1strain1:ko2strain1-ko1strain2:ko2strain2
PAIR = args.pair
# Fraction strain 1
FRACTION = args.fraction
# Unmodeled protein fraction
UNMODELED_PROTEIN = [args.unmodeled1, args.unmodeled2]
# Outer membrane transporter of crossfeed metabolite will be 65. * multiplier
SECRETION_KEFF_MULTIPLIER = [args.keff_transporter_1, args.keff_transporter_2]
MODE = args.mode
SCALE_SECRETION = str2bool(args.Scale_secretion)
GLUCOSE_UPTAKE = float(args.glucose_uptake)
DOCKER = str2bool(args.docker)

if MODE == 'default' or MODE == 'glucose_limited':
    RESTRICT = False
else:
    # RESTRICT can be 'experimental_inferred' or exchange reaction ID
    RESTRICT = args.Restrict_crossfeeding

print('PAIR', 'FRACTION', 'UNMODELED_PROTEIN', 'SCALE_SECRETION', 'RESTRICT')
print(PAIR, FRACTION, UNMODELED_PROTEIN, SCALE_SECRETION, RESTRICT)

# ************************************************************
# Load and prepare model for simulations
here = dirname(abspath(__file__))
resource_dir = resources.__path__[0]
with open('%s/iJL1678b_community.pickle' % resource_dir, 'rb') as f:
    model = pickle.load(f)

model.reactions.EX_glc__D_e_S1.lower_bound = GLUCOSE_UPTAKE / 2.
model.reactions.EX_glc__D_e_S2.lower_bound = GLUCOSE_UPTAKE / 2.
print(model.reactions.EX_glc__D_e_S2.lower_bound)
m_model = cobra.io.load_json_model('%s/iJO1366.json' % resource_dir)
print('loaded model')

# ============================================================
tic = time.time()


def solve_model(me_nlp, hs):
    # Re-compile expressions with new keffs in S and solve
    return me_nlp.bisectmu(precision=MU_PREC, mumin=MU_MIN, mumax=MU_MAX,
                           basis=hs)


def get_metabolic_flux(model):
    out_dict = {}
    for ex in model.reactions.query('EX_'):
        out_dict[ex.id] = ex.x
    for s in model.stoichiometric_data:
        for r in model.reactions.query(s.id + '_FWD_'):
            r_id = s.id + '_S1' if '_S1' in r.id else s.id + '_S2'
            out_dict[r_id] = r.x
        for r in model.reactions.query(s.id + '_REV_'):
            r_id = s.id + '_S1' if '_S1' in r.id else s.id + '_S2'
            out_dict[r_id] = -r.x
    return out_dict


kos1, kos2 = PAIR.split('-')

toc = time.time()-tic
print(kos1.split(':'), kos2.split(':'))

if MODE == 'secretion_keff_sweep':
    with open('%s/transport_reaction_per_strain.json' % resource_dir,
              'r') as f:
        transport_reaction = json.load(f)
    secretion_reactions = [transport_reaction[i] for i in PAIR.split('-')]
    secretion_keffs = [default_secretion_keff * i for i in SECRETION_KEFF_MULTIPLIER]
else:
    secretion_reactions = None
    secretion_keffs = None

setup_simulation(model, kos1.split(':'), kos2.split(':'), FRACTION,
                 unmodeled_protein_fractions=UNMODELED_PROTEIN,
                 secretion_reaction_keffs=secretion_keffs,
                 secretion_reactions=secretion_reactions,
                 restrict_uptake_flag=True,
                 scale_secretion=SCALE_SECRETION)

# Make sure strains cannot secrete any metabolites that it is auxotorphic for
for r in get_possible_uptake(m_model, kos1.split(':')):
    model.reactions.get_by_id(r + '_S1').upper_bound = 0
    print(str(kos1.split(':')), r)
for r in get_possible_uptake(m_model, kos2.split(':')):
    model.reactions.get_by_id(r + '_S2').upper_bound = 0
    print(str(kos2.split(':')), r)


# RESTRICT can either be "experimental_inferred" or exchange reaction ID
if RESTRICT:
    print('restricting exchange to one metabolite')
    with open('%s/restrict_exchange.json' % resource_dir, 'r') as f:
        restrict_exchange = json.load(f)

    for strain in [0, 1]:
        suffix = '_S1' if strain == 0 else '_S2'

        if RESTRICT == 'experimental_inferred':
            print("Using experimentally inferred exchange")
            restrict_reaction = restrict_exchange[[kos1, kos2][strain]][0]
        else:
            print('Restricting exchange to %s' % RESTRICT)
            restrict_reaction = ['EX_his__L_e', RESTRICT][strain]

        for r in model.reactions.query('EX_'):
            if not r.id.startswith('EX_') or '_Shared' in r.id or suffix not in r.id:
                continue
            if '_reverse' not in r.id:
                continue
            if model.reactions.get_by_id(
                    r.id.replace(suffix, '_Shared').replace(
                        '_reverse', '')).lower_bound < 0:
                continue
            if r.id != restrict_reaction + suffix + '_reverse':
                r.lower_bound = 0

    # Confirm proper reactions were limited
    for r in model.reactions.query('_reverse'):
        if r.lower_bound < 0:
            print(r.id, r.lower_bound, r.reaction)
    ITER = 0
else:
    ITER = 0

# ************************************************************
if MODE == 'unmodeled_sweep':
    out_directories = [os.getcwd(), MODE, PAIR,
                       '%.2f_%.2f_%i_unmodeled_protein' % (
                           UNMODELED_PROTEIN[0] * 100,
                           UNMODELED_PROTEIN[1] * 100,
                           float(ITER))]
elif MODE == 'secretion_keff_sweep':
    out_directories = [os.getcwd(), MODE, PAIR,
                       '%.2f_%.2f_%i_secretion_multiplier' % (
                           SECRETION_KEFF_MULTIPLIER[0] * 100,
                           SECRETION_KEFF_MULTIPLIER[1] * 100,
                           float(ITER))]
elif MODE == 'metabolite_limitation':
    out_directories = [os.getcwd(), MODE, PAIR, '%s' % RESTRICT]

elif MODE == 'default':
    out_directories = [os.getcwd(), MODE, PAIR,
                       '%.2f_%.2f_unmodeled_protein' % (
                           UNMODELED_PROTEIN[0] * 100,
                           UNMODELED_PROTEIN[1] * 100)]
elif MODE == 'glucose_limited':
    out_directories = [os.getcwd(), MODE, PAIR,
                       '%.2f_%.2f_unmodeled_protein' % (
                           UNMODELED_PROTEIN[0] * 100,
                           UNMODELED_PROTEIN[1] * 100)]

else:
    raise UserWarning('MODE not "unmodeled_sweep" or "secretion_keff_sweep" or'
                      ' "metabolite limitation" or "default" or "glucose_limited"')

output_dir = os.getcwd() if not DOCKER else '/output/'
for dir in out_directories[1:]:
    output_dir = os.path.join(output_dir, dir)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

# *************************************************************

with open('test_model.pickle', 'wb') as f:
    pickle.dump(model, f)

me_nlp = ME_NLP1(model, growth_key='mu')
me_nlp.compiled_expressions = me_nlp.compile_expressions()

# do not find community basis. This can cause issues
if os.path.isfile('community_basis_nope.pickle'):
    with open('community_basis.pickle', 'rb') as f:
        hs = pickle.load(f)
else:
    x, status, hs = me_nlp.solvelp(.0001)

muopt, hs, xopt, cache = solve_model(me_nlp, hs)

output_file = '%.2f_frac_strain1' % FRACTION
if model.solution is not None:
    model.solution.f = muopt
    with open(output_dir + '/' + output_file + '_sol.json', 'w') as f:
        json.dump(model.solution.x_dict, f)

    with open(output_dir + '/' + output_file + '_flux.json', 'w') as f:
        json.dump(get_metabolic_flux(model), f)

    with open('community_basis.pickle', 'wb') as f:
        pickle.dump(hs, f)

    rows = []
    for rxn_id, v in model.solution.x_dict.items():
        # Save rxn, flux, keff (if rxn has keff)
        rxn = model.reactions.get_by_id(rxn_id)
        keff = rxn.keff if hasattr(rxn, 'keff') else np.nan
        rows.append({'rxn': rxn_id, 'v': v, 'keff': keff})
    df_result = pd.DataFrame(rows)
else:
    df_result = pd.DataFrame([{'rxn': np.nan, 'v': np.nan,
                               'keff': np.nan}])
