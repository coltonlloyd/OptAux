import os
import tarfile

here = os.path.dirname(os.path.abspath(__file__))

for plot_kind in ['metabolite_limitation', 'secretion_keff_sweep', 'default']:
    tar_filename = '%s/community_me_sims/%s.tar.gz' % (here, plot_kind)
    if os.path.exists(tar_filename):
        tarfile.open(tar_filename).extractall('%s/community_me_sims' %
                                              here)