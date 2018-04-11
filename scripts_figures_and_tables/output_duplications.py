from __future__ import print_function, absolute_import, division

import os
import json

from optaux.ale_resequencing.find_duplications import (return_coverage_dict,
                                                       plot_coverage,
                                                       return_gene_duplicates)

here = os.path.dirname(os.path.abspath(__file__))
if __name__ == '__main__':
    alignment_loc = '/media/hard_drive/for_ale_analytics/aux/'
    save_loc = '%s/duplications' % here
    for pair in ['hisD_gltB', 'hisD_pyrC', 'hisD_gltA']:
        # coverage_dict = {'pair': {'ale': {'flask': {'isolate': {'replicate':
        # [coverage_per_position]}}}}}
        if not os.path.isfile('./%s_coverage_dict.json' % pair):
            coverage_dict = return_coverage_dict(alignment_loc, pair)
        else:
            print('using existing coverage dict')
            with open('%s_coverage_dict.json' % pair, 'r') as f:
                coverage_dict = json.load(f)

        duplicated_genes = return_gene_duplicates(pair, coverage_dict)

        plot_coverage(pair, coverage_dict, sections=50000, unfiltered=True)
