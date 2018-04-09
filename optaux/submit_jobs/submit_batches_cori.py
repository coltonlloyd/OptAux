import os
import numpy as np
from itertools import combinations
from optaux.resources import possible_uptake

#ko1s = [['HISTD'], ['HISTD'], ['HISTD']]

#ko2s = [['GLUDy', 'GLUSy'], ['DHORTS'], ['CS']]

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

modes = ['default']

submit_template = \
    "sbatch cori_submit_job %s %s %s %s %s --Restrict_crossfeeding %s " \
    "--keff_transporter_1 %s --keff_transporter_2 %s"

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
                        os.system(submit_template %
                                  (pair, f, q1, q2, mode,
                                   'experimental_inferred', k1, k2))

        elif mode == 'metabolite_limitation':
            for ko1, ko2 in zip(ko1s, ko2s):
                pair = ':'.join(ko1) + '-' + ':'.join(ko2)

                for met in possible_uptake.ko_uptakes[str(ko2)]:
                    q1 = .75
                    q2 = .75
                    k1 = 1.
                    k2 = 1.

                    for f in np.linspace(.05, .95, 10):
                        os.system(submit_template %
                                  (pair, f, q1, q2, mode, met, k1, k2))
        elif mode == 'unmodeled_sweep':
            for ko1, ko2 in zip(ko1s, ko2s):
                pair = ':'.join(ko1) + '-' + ':'.join(ko2)

                for q1, q2 in zip(np.linspace(.7, .8, 9),
                                  np.linspace(.8, .7, 9)):
                    k1 = 1.
                    k2 = 1.
                    for f in np.linspace(.05, .95, 10):
                        os.system(submit_template %
                                  (pair, f, q1, q2, mode,
                                   'experimental_inferred', k1, k2))

        elif mode == 'default':
            for ko1, ko2 in zip(ko1s, ko2s):
                pair = ':'.join(ko1) + '-' + ':'.join(ko2)
                q1 = .75
                q2 = .75
                k1 = 1.
                k2 = 1.
                for f in np.linspace(.05, .95, 10):
                    os.system(submit_template %
                              (pair, f, q1, q2, mode, False, k1, k2))
