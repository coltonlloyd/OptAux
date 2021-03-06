from __future__ import print_function, absolute_import

from os.path import abspath, dirname, exists

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from optaux.ale_resequencing import relative_abundance_plotting, \
    relative_abundance
from optaux.resources.update_mpl_rcparams import update_rcparams

update_rcparams()
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['legend.fontsize'] = 15

here = dirname(abspath(__file__))
if __name__ == '__main__':
    table_save_loc = '%s/relative_abundance/tables' % here
    fig_save_loc = '%s/relative_abundance/figures' % here
    alignment_loc = '/media/hard_drive/all_aux_ale_resequencing/aux/'  # alignment files not yet available

    fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharey=True, sharex='col')
    print(axes)
    for i, pair in enumerate(['hisDgltA', 'hisDgltB', 'hisDpyrC']):
        s1, s2 = pair[:4], pair[4:]
        char_file = \
            '%s/abundance_by_characteristic_%s.csv' % (table_save_loc, pair)
        if not exists(char_file):
            relative_abundance.get_characteristic_abundance_df(table_save_loc)
        char_df = pd.read_csv(char_file, index_col=0, header=None)

        cov_file = \
            '%s/abundance_by_coverage_%s.csv' % (table_save_loc, pair)
        if not exists(cov_file):
            relative_abundance.get_coverage_abundance_df(table_save_loc,
                                                         alignment_loc)
        cov_df = pd.read_csv(cov_file, index_col=0, header=None)

        relative_abundance_plotting.plot_abundance_df(char_df, s1, s2,
                                                      'characteristic',
                                                      axes[0][i])
        relative_abundance_plotting.plot_abundance_df(cov_df, s1, s2,
                                                      'coverage', axes[1][i])
    for i in [0, 1, 2]:
        axes[1][i].set_xlabel('ALE Number')

    for i in [0, 1]:
        axes[i][0].set_ylabel('Fractional Stain Abundance')

    fig.savefig('%s/coverage.eps' % fig_save_loc)
    fig.savefig('%s/coverage.png' % fig_save_loc)

    # Normalize = True ensure abundance predictions add to 1
    # NOTE: Filter for unique case where coverage is too low to identify
    # characteristic mutation, so filter out ale_10_flask_23
    relative_abundance_plotting.plot_pairwise_comparison_of_preditions(
        table_save_loc, fig_save_loc, normalize=True,
        filter_ale_10_flask_23=True)
