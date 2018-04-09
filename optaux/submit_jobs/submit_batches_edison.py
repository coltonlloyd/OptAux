import os
ko1s = [['CS'], ['CS'], ['CS'], ['HISTD'], ['HISTD'], ['DHORTS']]

ko2s = [['GLUDy', 'GLUSy'], ['HISTD'], ['DHORTS'],
        ['GLUDy', 'GLUSy'], ['DHORTS'], ['GLUDy', 'GLUSy']]
for ko1, ko2 in zip(ko1s, ko2s):
    pair = ':'.join(ko1) + '-' + ':'.join(ko2)
    for q in [.7, .8, .9, 1., 1.1, 1.2, 1.3]:
        source_dir = os.getcwd() + '/community_sims'
        if not os.path.isdir(source_dir):
            os.mkdir(source_dir)
        if not os.path.isdir(os.path.join(source_dir, pair)):
            os.mkdir(os.path.join(source_dir, pair))
        unmodled_name = '%i_unmodeled_protein' % (q * 100.)
        if not os.path.isdir(os.path.join(source_dir, pair, unmodled_name)):
            os.mkdir(os.path.join(source_dir, pair, unmodled_name))

        os.system("sbatch edison_submit_job %s --Unmodeled_protein %s" % (pair,
                                                                          q))
