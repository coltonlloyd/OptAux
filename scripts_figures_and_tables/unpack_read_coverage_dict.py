import os
import tarfile

here = os.path.dirname(os.path.abspath(__file__))

for pair in ['hisD_gltA', 'hisD_gltB', 'hisD_pyrC']:
    tar_filename = '%s/duplications/%s_coverage_dict.tar.xz' % (here, pair)
    if os.path.exists(tar_filename):
        tarfile.open(tar_filename).extractall('%s/duplications' %
                                              here)