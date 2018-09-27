import os
import tarfile

here = os.path.dirname(os.path.abspath(__file__))

for model in ['null', 'default']:
    for plot_kind in ['metabolite_limitation', 'secretion_keff_sweep',
                      'default', 'glucose_limited']:
        tar_filename = '%s/community_sims_output_%s_keffs/%s.tar.gz' % (here, model, plot_kind)
        if os.path.exists(tar_filename):
            tarfile.open(tar_filename).extractall('%s/community_sims_output_%s_keffs' %
                                                  (here, model))