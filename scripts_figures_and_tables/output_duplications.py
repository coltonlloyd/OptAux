from __future__ import print_function, absolute_import, division

import os
import json

from optaux.ale_resequencing.find_duplications import (return_coverage_dict,
                                                       plot_coverage,
                                                       return_gene_duplicates)

here = os.path.dirname(os.path.abspath(__file__))
if __name__ == '__main__':
    alignment_loc = ''  # alignment files not yet available
    save_loc = '%s/duplications/' % here
    for pair in ['hisD_gltB', 'hisD_pyrC', 'hisD_gltA']:
        # coverage_dict = {'pair': {'ale': {'flask': {'isolate': {'replicate':
        # [coverage_per_position]}}}}}
        if os.path.isfile('%s/%s_coverage_dict.json' % (save_loc, pair)):
            print('using existing coverage dict')
            with open('%s/%s_coverage_dict.json' % (save_loc, pair), 'r') as f:
                coverage_dict = json.load(f)

        else:
            coverage_dict = return_coverage_dict(save_loc, alignment_loc, pair)

        duplicated_genes = \
            return_gene_duplicates('%s/duplicated_genes' % save_loc, pair,
                                   coverage_dict)


        plot_coverage(save_loc, pair, coverage_dict, sections=50000,
                      unfiltered=True)
