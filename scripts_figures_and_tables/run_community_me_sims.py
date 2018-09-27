from __future__ import print_function, division, absolute_import

import os

import numpy as np

from optaux.submit_jobs.submit_local import does_sim_exist, run_pool, \
    submit_job
from optaux.resources import possible_uptake

ko1s = [['HISTD'], ['HISTD'], ['HISTD']]

ko2s = [['CS'], ['GLUDy', 'GLUSy'], ['DHORTS']]

modes = ['secretion_keff_sweep', 'metabolite_limitation',
         'default', 'glucose_limited']


here = os.path.dirname(os.path.abspath(__file__))
simulation_directory = '%s/community_sims_output_null_keffs' % here

# The communities cannot grow when exchanging these metabolites
skip_mets = ['EX_fe3dcit_e', 'EX_fe3dcit_e', 'EX_progly_e', 'EX_23ccmp_e',
             'EX_23cump_e', 'EX_3cmp_e', 'EX_3ump_e', 'EX_udpacgal_e']

if __name__ == '__main__':
    values = []
    for model in ['default', 'null']:
        for mode in modes:
            if mode == 'secretion_keff_sweep':
                for ko1, ko2 in zip(ko1s, ko2s):
                    pair = ':'.join(ko1) + '-' + ':'.join(ko2)

                    for k1, k2 in zip(np.linspace(-2, 2, 9),
                                      np.linspace(2, -2, 9)):
                        q1 = q2 = .6
                        k1 = round(10**k1, 4)
                        k2 = round(10**k2, 4)

                        for f in np.linspace(.05, .95, 10):
                            if not does_sim_exist(simulation_directory, pair, mode,
                                                  f, q1, q2, k1, k2):
                                values.append((pair, f, q1, q2, mode,
                                               'experimental_inferred', k1, k2,
                                               -1000, model))

            elif mode == 'metabolite_limitation':
                for ko1, ko2 in zip(ko1s, ko2s):
                    pair = ':'.join(ko1) + '-' + ':'.join(ko2)

                    for met in possible_uptake.ko_uptakes[str(ko2)]:
                        q1 = q2 = .6
                        k1 = k2 = 1.

                        # Skip if community cannot grow with the metabolite
                        if met in skip_mets:
                            continue

                        for f in np.linspace(.05, .95, 10):
                            if not does_sim_exist(simulation_directory, pair, mode,
                                                  f, q1, q2, k1, k2, met):
                                values.append(
                                    (pair, f, q1, q2, mode, met,  k1, k2, -1000,
                                     model))
            elif mode == 'unmodeled_sweep':
                for ko1, ko2 in zip(ko1s, ko2s):
                    pair = ':'.join(ko1) + '-' + ':'.join(ko2)

                    for q1, q2 in zip(np.linspace(.7, .8, 9),
                                      np.linspace(.8, .7, 9)):
                        k1 = k2 = 1.

                        for f in np.linspace(.05, .95, 10):
                            if not does_sim_exist(simulation_directory, pair, mode,
                                                  f, q1, q2, k1, k2):
                                values.append((pair, f, q1, q2, mode,
                                               'experimental_inferred', k1, k2,
                                               -1000, model))

            elif mode == 'default':
                for ko1, ko2 in zip(ko1s, ko2s):
                    pair = ':'.join(ko1) + '-' + ':'.join(ko2)
                    q1 = q2 = .6
                    k1 = k2 = 1.
                    for f in np.linspace(.05, .95, 10):
                        if not does_sim_exist(simulation_directory, pair, mode, f,
                                              q1, q2, k1, k2):
                            values.append((pair, f, q1, q2, mode, False, k1, k2,
                                           -1000, model))

            elif mode == 'glucose_limited':
                for ko1, ko2 in zip(ko1s, ko2s):
                    pair = ':'.join(ko1) + '-' + ':'.join(ko2)
                    q1 = .36
                    q2 = .36
                    k1 = 1.
                    k2 = 1.
                    for f in np.linspace(.05, .95, 10):
                        if not does_sim_exist(simulation_directory, pair, mode, f,
                                              q1, q2, k1, k2, model):
                            values.append((pair, f, q1, q2, mode, False, k1, k2,
                                           -5, model))
        print(values)
        print(len(values))
    #run_pool(submit_job, values, processes=4)
