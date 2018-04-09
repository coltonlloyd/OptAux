import os
import numpy as np
from torque_create_job import apply_torque
ko1s = [['CS'], ['CS'], ['CS'], ['HISTD'], ['HISTD'], ['DHORTS']]

ko2s = [['GLUDy', 'GLUSy'], ['HISTD'], ['DHORTS'],
        ['GLUDy', 'GLUSy'], ['DHORTS'], ['GLUDy', 'GLUSy']]

ko1s = [['HISTD']]
ko2s = [['DHORTS']]


source_dir = os.getcwd() + '/community_sims'
if not os.path.isdir(source_dir):
    os.mkdir(source_dir)

for ko1, ko2 in zip(ko1s, ko2s):
    pair = ':'.join(ko1) + '-' + ':'.join(ko2)

    if not os.path.isdir(os.path.join(source_dir, pair)):
        os.mkdir(os.path.join(source_dir, pair))
    for q in [.4]:#, .8, .9, 1., 1.1, 1.2, 1.3]:
        unmodeled_name = '%i_%i_%i_unmodeled_protein' % (q * 100., q * 100., 0)
        if not os.path.isdir(os.path.join(source_dir, pair, unmodeled_name)):
            os.mkdir(os.path.join(source_dir, pair, unmodeled_name))

        for f in np.linspace(.1, .9, 10):
            apply_torque(pair, f, q, q)