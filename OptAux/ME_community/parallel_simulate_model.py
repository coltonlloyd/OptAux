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
from mpi4py import MPI
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
from OptAux.ME_community.run_ALE_pairs import setup_simulation

# ************************************************************
# Parameters
# ------------------------------------------------------------
# Bisection parameters
MU_PREC = 1e-3
MU_MIN = 0.
MU_MAX = 2.5

simstr = ''

# ************************************************************
# Argument
parser = argparse.ArgumentParser(description='Simulation parameters.')

parser.add_argument('pair', help='', type=str)
parser.add_argument('--N_sims', help='', type=int, default=10)
parser.add_argument('--Start',  help='', default=.1)
parser.add_argument('--Stop', help='', default=.9)
parser.add_argument('--Scale_exchange', default=False)
parser.add_argument('--Restrict_crossfeeding', default=False)
parser.add_argument('--Unmodeled_protein', default=1.)

args = parser.parse_args()

# pair must be in form of ko1strain1:ko2strain1-ko1strain2:ko2strain2
PAIR = args.pair
N_SIMULATIONS = int(args.N_sims)
START_VALUE = float(args.Start)
STOP_VALUE = float(args.Stop)
SCALE = args.Scale_exchange
RESTRICT = args.Restrict_crossfeeding
UNMODELED_PROTEIN = float(args.Unmodeled_protein)

# ************************************************************
# Load and prepare model for simulations
here = dirname(abspath(__file__))
resources = relpath('../resources', here)
with open('%s/iJL1678b_community.pickle' % resources, 'r') as f:
    model = pickle.load(f)

print('loaded model')
# ************************************************************

# Get MPI worker information
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# number of cores running process
nWorkers = size
# ------------------------------------------------------------

output_dir = PAIR
if not os.path.isdir(os.path.join(os.getcwd(), 'community_sims', output_dir)):
    os.mkdir(os.path.join(os.getcwd(), 'community_sims', output_dir))


# ------------------------------------------------------------
# Dictionary of work:
samples = range(0, N_SIMULATIONS)

# Used for parameter sweeps / growth curves
value_list = np.linspace(START_VALUE, STOP_VALUE, N_SIMULATIONS)

work_dict_list = [{'sample': sample} for sample in samples]
inds_all = range(0, len(work_dict_list))
# -----------------------------------------------------------

###
print(nWorkers, 'workers available')
print('%d tasks' % len(inds_all))
###

# ------------------------------------------------------------
# Split the work
nTasks = len(inds_all)
tasks_per_worker = nTasks/nWorkers
rem = nTasks - tasks_per_worker * nWorkers
# Distributes remainder across workers
worker_sizes = np.array([tasks_per_worker + (1 if i < rem else 0)
                         for i in range(0, nWorkers)], 'i')

# ------------------------------------------------------------
# Thus, worker_tasks[rank] gives the actual indices to work on
inds_end = np.cumsum(worker_sizes)
inds_start = np.hstack((0, inds_end[0:-1]))
worker_tasks = [inds_all[inds_start[i]:inds_end[i]] for i in range(0,
                                                                   nWorkers)]

print('inds_all:', inds_all)
print('worker_tasks:', worker_tasks)
print('%d batches' % len(worker_tasks))

# ============================================================
# Work performed by each worker

# Compile Expressions
me_nlp = ME_NLP1(model, growth_key='mu')
me_nlp.compiled_expressions = me_nlp.compile_expressions()


def do_work(inds):
    # Solve
    def simulate(ind, hs):
        # ----------------------------------------------------
        # Choose sample
        sim_dict = work_dict_list[ind]
        sample = sim_dict['sample']
        print('Simulation sample %d in job %s' % (sample, rank))
        # ----------------------------------------------------
        tic = time.time()

        def solve_model(model):
            # Re-compile expressions with new keffs in S and solve
            me_nlp = ME_NLP1(model, growth_key='mu')
            return me_nlp.bisectmu(precision=MU_PREC, mumin=MU_MIN, mumax=MU_MAX,
                                   basis=hs)

        value = value_list[sample]
        kos1, kos2 = PAIR.split('-')

        toc = time.time()-tic
        setup_simulation(model, kos1.split(':'), kos2.split(':'), value,
                         UNMODELED_PROTEIN, restrict_uptake_flag=RESTRICT,
                         scale_secretion_and_uptake=SCALE)

        muopt, hs, xopt, cache = solve_model(model)

        output_file = '%.2f_frac_strain1.pickle' % value
        if model.solution is not None:
            with open(output_dir + '/' + output_file + '_flux.json', 'wb') as f:
                json.dump(model.get_metabolic_flux(), f)

            with open(output_dir + '/' + output_file + '_sol.pickle', 'wb') as f:
                pickle.dump(model.solution, f)

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

        # ----------------------------------------------------
        # Write result to file
        #timestr = time.strftime("%Y%m%d")
        #filename = 'result_sampleme%d%s_%s_job%s_%d.csv' % (
        #    ME_PROTO_VER, simstr,
        #    timestr, SIMULATION_TYPE, sample)

        # Save solution
        #df_result.to_csv(filename, index=False)

        # Return
        result = {'basis': hs, 'xopt': xopt, 'muopt': muopt}
        return result, df_result, toc

    # Use correct indices and allow warm-start within a worker if needed
    sols = []

    hs = None
    for ind in inds:
        sol = simulate(ind, hs)
        res = sol[0]
        hs_new = res['basis']
        sols.append(sol)
        if hs_new is not None:
            hs = hs_new

    df_results = [sol[1] for sol in sols]
    times = [sol[2] for sol in sols]

    results = {'ind': inds, 'time': times}

    return results

# ============================================================
# Do work on data chunk
data = do_work(worker_tasks[rank])
print('Finished work by worker %d' % rank)

# ------------------------------------------------------------
# Gather results by root
data = comm.gather(data, root=0)

# ------------------------------------------------------------
# Report final result
if rank == 0:
    # Gather results
    print('Gathered data by root')

    # Save to pickle first
    timestr = time.strftime("%Y%m%d_%H%M")
    filename = 'results_community_%.2f_unmodeled_%s_%s_job%s.pickle' % (
        UNMODELED_PROTEIN, simstr, timestr, PAIR)

    with open(filename, 'wb') as iofile:
        pickle.dump(data, iofile)

    print('Saved gathered results to pickle file:', filename)

    # Save as dataframe
    inds = [i for res in data for i in res['ind']]
    times = [t for res in data for t in res['time']]

    df = pd.DataFrame({'ind': inds, 'time': times})

    filename_df = 'results_community_%.2f_unmodeled_%s_job%s_%s.csv' % (
        UNMODELED_PROTEIN, simstr, timestr, PAIR)
    df.to_csv(filename_df, index=False)

    print('Saved results to csv:', filename_df)
