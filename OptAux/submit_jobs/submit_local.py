from multiprocessing import Pool
import os
from os.path import abspath, dirname
import numpy as np
from OptAux import ME_community
from OptAux.resources import possible_uptake
import subprocess
from itertools import combinations


ko1s = [['HISTD'], ['HISTD'], ['HISTD']]

ko2s = [['CS'], ['GLUDy', 'GLUSy'], ['DHORTS']]

# From Mee et al.
ko_to_rxns = {'hisB': ['HISTP', 'IGPDH'],
              'cysE': ['SERAT'],
              'pheA': ['PPNDH'],
              'lysA': ['DAPDC'],
              'leuB': ['IPMD'],
              'metA': ['HSST'],
              'argA': ['ACGS'],
              'thrC': ['THRS', '4HTHRS'],
              'trpC': ['PRAIi', 'IGPS'],
              'tyrA': ['PPND']}

ko2s = []
ko1s = []
for kos1, kos2 in combinations(list(ko_to_rxns.values()), 2):
    ko1s.append(kos1)
    ko2s.append(kos2)

ko2s = ko2s[0:1]
ko1s = ko1s[0:1]
#ko1s = [['HISTD']]
#ko2s = [['CS']]
#modes = ['unmodeled_sweep', 'secretion_keff_sweep', 'metabolite_limitation']
modes = ['default']

sim_script_dir = ME_community.__path__[0]
simulation_directories = '/home/sbrg-cjlloyd/mee_aux_compare'

# The communities cannot grow when exchanging these metabolites
skip_mets = ['EX_fe3dcit_e', 'EX_fe3dcit_e', 'EX_progly_e', 'EX_23ccmp_e',
             'EX_23cump_e', 'EX_3cmp_e', 'EX_3ump_e', 'EX_udpacgal_e']


def run_pool(function, values, processes=4):
    pool = Pool(processes=processes)
    pool.map(function, values)


def submit_job(value_tuple):
    pair, f, q1, q2, mode, restrict, k1, k2 = value_tuple
    os.system("python2.7 %s/simulate_model.py %s %s %s %s %s "
              "--Restrict_crossfeeding %s --keff_transporter_1 %s "
              "--keff_transporter_2 %s" %
              (sim_script_dir, pair, f, q1, q2, mode, restrict, k1, k2))


def does_sim_exist(pair, mode, f, q1, q2, k1, k2, met=None):
    if mode == 'secretion_keff_sweep':
        file_loc = \
            '%s/%s/%s/%.2f_%.2f_0_secretion_multiplier/%.2f_frac_strain1_sol.pickle' % (simulation_directories, mode, pair, k1*100, k2*100, f)

    elif mode == 'unmodeled_sweep':
        file_loc = \
            '%s/%s/%s/%.2f_%.2f_0_unmodeled_protein/%.2f_frac_strain1_sol.pickle' % (simulation_directories, mode, pair, q1*100, q2*100, f)

    elif mode == 'metabolite_limitation':
        file_loc = \
            '%s/%s/%s/%s/%.2f_frac_strain1_sol.pickle' % \
            (simulation_directories, mode, pair, met, f)

    elif mode == 'default':
        file_loc = \
            '%s/%s/%s/%.2f_%.2f_unmodeled_protein/%.2f_frac_strain1_sol.pickle' % (simulation_directories, mode, pair, q1*100, q2*100, f)

    if os.path.exists(file_loc):
        return True
    else:
        print(file_loc)
        return False


if __name__ == '__main__':
    values = []
    for mode in modes:
        if mode == 'secretion_keff_sweep':
            for ko1, ko2 in zip(ko1s, ko2s):
                pair = ':'.join(ko1) + '-' + ':'.join(ko2)

                for k1, k2 in zip(np.linspace(-2, 2, 9),
                                  np.linspace(2, -2, 9)):
                    q1 = .75
                    q2 = .75
                    k1 = round(10**k1, 4)
                    k2 = round(10**k2, 4)

                    for f in np.linspace(.05, .95, 10):
                        if not does_sim_exist(pair, mode, f, q1, q2, k1, k2):
                            values.append((pair, f, q1, q2, mode,
                                           'experimental_inferred', k1, k2))

        elif mode == 'metabolite_limitation':
            for ko1, ko2 in zip(ko1s, ko2s):
                pair = ':'.join(ko1) + '-' + ':'.join(ko2)

                for met in possible_uptake.ko_uptakes[str(ko2)]:
                    q1 = .75
                    q2 = .75
                    k1 = 1.
                    k2 = 1.

                    # Skip if community cannot grow with the metabolite
                    if met in skip_mets:
                        continue

                    for f in np.linspace(.05, .95, 10):
                        if not does_sim_exist(pair, mode, f, q1, q2, k1, k2,
                                              met):
                            values.append(
                                (pair, f, q1, q2, mode, met,  k1, k2))
        elif mode == 'unmodeled_sweep':
            for ko1, ko2 in zip(ko1s, ko2s):
                pair = ':'.join(ko1) + '-' + ':'.join(ko2)

                for q1, q2 in zip(np.linspace(.7, .8, 9),
                                  np.linspace(.8, .7, 9)):
                    k1 = 1.
                    k2 = 1.
                    for f in np.linspace(.05, .95, 10):
                        if not does_sim_exist(pair, mode, f, q1, q2, k1, k2):
                            values.append((pair, f, q1, q2, mode,
                                           'experimental_inferred', k1, k2))

        elif mode == 'default':
            for ko1, ko2 in zip(ko1s, ko2s):
                pair = ':'.join(ko1) + '-' + ':'.join(ko2)
                q1 = .75
                q2 = .75
                k1 = 1.
                k2 = 1.
                for f in np.linspace(.05, .95, 1):
                    if not does_sim_exist(pair, mode, f, q1, q2, k1, k2):
                        values.append((pair, f, q1, q2, mode,
                                       False, k1, k2))

    run_pool(submit_job, values, processes=2)
