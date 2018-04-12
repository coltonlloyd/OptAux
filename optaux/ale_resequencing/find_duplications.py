from __future__ import print_function

import matplotlib
matplotlib.use('Agg')
import json
from Bio import SeqIO
from collections import Counter, defaultdict
import pysam
import cobra.test
import seaborn as sns
import os
import glob
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

from optaux import resources

plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
plt.rcParams['axes.labelsize'] = 35
plt.rcParams['axes.titlesize'] = 30
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['grid.color'] = 'white'
plt.rcParams['figure.edgecolor'] = 'black'
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['legend.fontsize'] = 15.0
plt.rcParams['axes.grid'] = False
plt.rcParams['font.size'] = 20

iJO1366 = cobra.test.create_test_model('ecoli')

resource_dir = resources.__path__[0]
genbank = '%s/resequencing_data/CP009273_1.gb' % resource_dir
with open('%s/gene_name_to_gb_info.json' % resource_dir, 'r') as f:
    g_to_info = json.load(f)


def get_gene_position_info(genbank_loc):
    seqs = SeqIO.read(genbank_loc, 'gb')
    a = {}
    gene_to_pos = {}
    for feat in seqs.features:
        if feat.type != 'CDS':
            continue
        a[feat.qualifiers['gene'][0]] = '%i-%i' % (
        feat.location.start, feat.location.end)
    for key, value in a.items():
        start, end = value.split('-')
        gene_to_pos[tuple(range(int(start), int(end)))] = key

    return gene_to_pos


def get_average_and_stdev_of_sections(coverage_list, genome_size,
                                      num_sections=10000):
    """Divides genomes into <num_sections> fragments and returns the average
    and standard deviations of coverages for each of the fragments, rounded to
    nearest 10"""

    i_start = 0
    mean_list = []
    std_list = []
    for i, stop in enumerate(np.linspace(1, genome_size, num_sections)):
        if i == 0:
            continue
        section = np.array(coverage_list[i_start: int(stop)])
        mean = section.mean()
        std = section.std()
        if not np.isnan(mean):
            mean_list.append(round(mean, -1))
        if not np.isnan(std):
            std_list.append(round(std, -1))
        i_start = int(stop)
    return mean_list, std_list


def return_summary_file_loc(ale, flask, isolate, replicate):
    path = '%s/resequencing_data/%s-%s-%s-%s/summary.json' % \
           (resource_dir, ale, flask, isolate, replicate)
    if os.path.exists(path):
        return path
    else:
        raise UserWarning(path, 'not valid')


def get_fit_mean_from_breseq(path):
    print(path)
    with open(path, 'r') as f:
        summary = json.load(f)
    # The replicate is either 1 or 2 try both

    fit_mean = float(summary['coverage_nbinom_mean_parameter'])
    return fit_mean


def find_duplicated_genes(duplicated_genes, coverage_list, mean, pair, ale_num,
                          flask_num, isolate, replicate, cutoff):
    """
    If > 80% of gene basepairs are above the 1.25 * mean cutoff, then this
    gene is considered duplicated
    """

    df = pd.DataFrame()

    gene_to_pos = get_gene_position_info(genbank)

    # Determine which genes duplicate in a group
    group = 0
    last = 0
    for index, pos in enumerate(sorted(gene_to_pos)):
        gene = gene_to_pos[pos]
        pos = set(pos)

        # Make list of positions within gene with coverage above threshold
        dup_pos = []
        for p in pos:
            if coverage_list[p] > (cutoff * mean):
                dup_pos.append(p)

        # Append gene if >80% contains duplicated base pair
        if len(pos.intersection(dup_pos)) > (.8 * len(pos)):
            multiplicity = []
            # Iterate to new group if genes are not located next to each other
            if index - last > 5:
                group += 1
            duplicated_genes[ale_num][flask_num].append(gene)
            for i in pos:
                multiplicity.append(coverage_list[i])
            df.loc[gene, 'Multiplicity'] = np.array(multiplicity).mean() / mean
            df.loc[gene, 'Group Num'] = group
            df.loc[gene, 'Start Position'] = list(pos)[0]
            try:
                df.loc[gene, 'Product'] = g_to_info[gene]['gene_product']
            except:
                df.loc[gene, 'Product'] = ''
            try:
                df.loc[gene, 'Gene Name'] = g_to_info[gene]['gene_name']
            except:
                df.loc[gene, 'Gene Name'] = ''
            try:
                df.loc[gene, 'Reactions'] = str(iJO1366.genes.get_by_id(df.loc[gene, 'Gene Name']).reactions)
            except:
                df.loc[gene, 'Reactions'] = ''
            last = index
    df.to_csv('./dups/%s_%s_%s_%s_%s_genes.csv' %
              (pair, ale_num, flask_num, isolate, replicate))
    return duplicated_genes


def filter_starting_strain_dups(coverage_dict, pair, ale):

    starting_dup_pos = np.zeros(len(coverage_dict[ale]["0"]["0"]["1"]))
    for isolate in ["0", "1"]:
        print(len(coverage_dict[ale]["0"][isolate]['1']))
        starting_coverage = coverage_dict[ale]["0"][isolate]['1']
        file = return_summary_file_loc(pair, ale, "0", int(isolate), "1")
        mean = get_fit_mean_from_breseq(file)
        array = np.array([1 if i > (2 * mean) else 0 for i in
                          starting_coverage])

        starting_dup_pos += array

    return starting_dup_pos


def return_gene_duplicates(pair, coverage_dict, cutoff=1.25):
    all_dup_genes = {}

    all_dup_genes[pair] = {}
    starting_dup_pos = \
        filter_starting_strain_dups(coverage_dict, pair, "0")

    for ale in coverage_dict:
        print(pair, ale)
        all_dup_genes[ale] = defaultdict(list)
        for flask in coverage_dict[ale]:
            for isolate in coverage_dict[ale][flask]:
                if isolate not in coverage_dict[ale][flask]:
                    print(ale, flask, isolate)
                    continue
                for replicate, read_coverage in coverage_dict[ale][flask][isolate].items():
                    print(pair, ale, flask, isolate, replicate)
                    summary_file = return_summary_file_loc(ale, flask,
                                                           isolate, replicate)
                    mean = get_fit_mean_from_breseq(summary_file)

                    read_cov = []
                    for i, cov in enumerate(read_coverage):
                        filtered_cov = cov if starting_dup_pos[i] == 0 else mean
                        read_cov.append(filtered_cov)

                    find_duplicated_genes(all_dup_genes, read_cov, mean,
                                          pair, ale, flask, isolate, replicate,
                                          cutoff)
    return all_dup_genes


def return_coverage_dict(alignment_loc, pair, genome_size=4631468):
    """Return dictionary of {ale:{flask:{isolate:{replicate: coverage_list}}}

    where coverage_list is a list of integers corresponding to the number of
    alignments at each base pair postition

    """
    all_coverage = {}

    for file in glob.glob(alignment_loc +
                          'aux_%s/breseq/ale/*/data/reference.bam' % pair):

        # Get the ale information
        ale_name = [str(int(i)) for i in
                    file.split('ale/')[1].split('/data')[0].split('-')]
        ale_num, flask, isolate, replicate = ale_name

        def _add_ale_info_to_coverage_dict(ale_num, flask, isolate,
                                           all_coverage):
            if ale_num not in all_coverage:
                all_coverage[ale_num] = {}

            if flask not in all_coverage[ale_num]:
                all_coverage[ale_num][flask] = {}

            if isolate not in all_coverage[ale_num][flask]:
                all_coverage[ale_num][flask][isolate] = {}

        _add_ale_info_to_coverage_dict(ale_num, flask, isolate, all_coverage)

        # Get coverage for each basepair
        samfile = pysam.AlignmentFile(file, 'rb')
        read_cov_dict = {}
        for pileupcolumn in samfile.pileup('CP009273', int(0),
                                           int(genome_size)):
            read_cov_dict[pileupcolumn.pos] = pileupcolumn.n
        samfile.close()

        # If coverage at a position is 0, there is no pileup() entry
        read_cov = []
        start = -1
        for key in sorted(read_cov_dict):
            diff = key - start
            if diff == 1:
                read_cov.append(read_cov_dict[key])
            else:
                # zero instead of one
                while diff > 0:
                    read_cov.append(0)
                    diff -= 1
            start = key
        print(len(read_cov))
        print(pair, ale_num, flask, isolate)
        all_coverage[ale_num][flask][isolate][replicate] = read_cov

    with open('%s_coverage_dict.json' % pair, 'w') as f:
        json.dump(all_coverage, f)

    return all_coverage


def add_gene_positions(ax, ymax, kind='all'):
    if kind == 'all':
        pass

        #
    # Prophages
    elif kind == 'endclone':
        ax.plot([258669, 258669], [0, ymax], 'k--', linewidth=.5)  # CP4-6
        ax.plot([292976, 292976], [0, ymax], 'k--', linewidth=.5)  # CP4-6

        ax.plot([560258, 560258], [0, ymax], 'k--', linewidth=.5)  # DLP12
        ax.plot([581559, 581559], [0, ymax], 'k--', linewidth=.5)  # DLP12

        ax.plot([1191676, 1191676], [0, ymax], 'k--', linewidth=.5)  # e14
        ax.plot([1206868, 1206868], [0, ymax], 'k--', linewidth=.5)  # e14

        ax.plot([1406156, 1406156], [0, ymax], 'k--', linewidth=.5)  # Rac
        ax.plot([1429215, 1429215], [0, ymax], 'k--', linewidth=.5)  # Rac

        ax.plot([1626543, 1626543], [0, ymax], 'k--', linewidth=.5)  # Qin/kim
        ax.plot([1647000, 1647000], [0, ymax], 'k--', linewidth=.5)  # Qin/kim

        ax.plot([2060835, 2060835], [0, ymax], 'k--', linewidth=.5)  # cp4-44
        ax.plot([2072512, 2072512], [0, ymax], 'k--', linewidth=.5)  # cp4-44

        ax.plot([2160679, 2160679], [0, ymax], 'k--', linewidth=.5)  # pr-x
        ax.plot([2161307, 2161307], [0, ymax], 'k--', linewidth=.5)  # pr-x

        ax.plot([2459864, 2459864], [0, ymax], 'k--', linewidth=.5)  # cps-53
        ax.plot([2470078, 2470078], [0, ymax], 'k--', linewidth=.5)  # cps-53

        ax.plot([2552058, 2552058], [0, ymax], 'k--', linewidth=.5)  # Eut/CPZ-55
        ax.plot([2558820, 2558820], [0, ymax], 'k--', linewidth=.5)  # Eut/CPZ-55

        ax.plot([2749316, 2749316], [0, ymax], 'k--', linewidth=.5)  # CP4-57
        ax.plot([2771345, 2771345], [0, ymax], 'k--', linewidth=.5)  # CP4-57

    elif kind == 'gltJ':
        ax.text(1394503, 24.5, 'abgT\n$\downarrow$', ha='center',
                va='top', fontdict={'fontsize': 20})
        ax.text(681984, 24.5, 'gltJ\n$\downarrow$', ha='center',
                va='top', fontdict={'fontsize': 20})
    elif kind == 'hisJ':
        ax.text(2419484, 4.8, 'hisJ\n$\downarrow$', ha='center',
                va='top', fontdict={'fontsize': 20})


def plot_and_format(ax, x, y, genome_size, ale, flask, isolate, replicate,
                    mean, xrange=None):

    xrange = [0, genome_size] if not xrange else xrange

    ymax = max(y / mean) if max(y / mean) > 2 else 2

    ax.plot(x, y / mean, linewidth=2.5)

    ax.set_ylim([0, ymax])
    ax.set_xlim(xrange)

    ax.set_title('A%i F%i I%i R%i' %
                 (ale, flask, isolate, replicate))
    ax.locator_params(axis='y', nbins=4)
    ax.legend(loc='upper center')
    ax.set_facecolor('#E0E0E0')


def save_fig(fig, pair, suffix, unfiltered=False):
    fig.text(0.5, -0.02, 'Genome position', ha='center', size=35)
    fig.text(-0.01, 0.5, 'Multiplicity', va='center', rotation='vertical',
             size=35)
    fig.tight_layout()
    if unfiltered:
        fig.savefig(os.path.expanduser('~/Dropbox/unfiltered_%s_%s.png' % (pair, suffix)),
                    bbox_inches='tight')
        fig.savefig(os.path.expanduser('~/Dropbox/unfiltered_%s_%s.svg' % (pair, suffix)),
                    bbox_inches='tight')
    else:
        fig.savefig(os.path.expanduser('~/Dropbox/filtered_%s_%s.png' % (pair, suffix)),
                    bbox_inches='tight')
        fig.savefig(os.path.expanduser('~/Dropbox/filtered_%s_%s.svg' % (pair, suffix)),
                    bbox_inches='tight')


def plot_coverage(pair, cov_dict, sections=10000, cutoff=1.25,
                  genome_size=4631468, unfiltered=False):

    # change string keys to integer keys
    coverage_dict = {}

    # There are a variable number of resequencing data for each ALE. Keep
    # track of the maximum number of experiments
    max_cols = 0
    for ale in cov_dict:
        coverage_dict[int(ale)] = {}
        cols = 0
        for flask in cov_dict[ale]:
            coverage_dict[int(ale)][int(flask)] = {}
            for isolate in cov_dict[ale][flask]:
                coverage_dict[int(ale)][int(flask)][int(isolate)] = {}
                for replicate, cov in cov_dict[ale][flask][isolate].items():
                    cols += 1
                    if cols > max_cols:
                        max_cols = cols
                    coverage_dict[int(ale)][int(
                        flask)][int(isolate)][int(replicate)] = cov

    plots = 2
    fig_gltJ, axes_gltJ = plt.subplots(plots, max_cols, figsize=(4*max_cols,
                                                                 3 * plots))

    # To clean up output plots, remove duplications seen in starting strains.
    # These correspond to IS elements, etc.
    starting_dup_pos = \
        filter_starting_strain_dups(cov_dict, pair, "0")

    print(pair, 'max cols ', max_cols)
    plots = 4 if pair == 'hisD_gltA' else 3
    fig, axes = plt.subplots(max_cols, plots,
                             figsize=(7 * plots, 3 * max_cols),
                             sharex=True)
    plt.locator_params(nbins=4)

    a = -1
    for ale in sorted(coverage_dict):
        if ale in [7, 1, 0]:
            continue
        a += 1
        print(pair, ale)
        num_cols = -1
        for flask in sorted(coverage_dict[ale]):
            for isolate in sorted(coverage_dict[ale][flask]):
                for replicate in sorted(coverage_dict[ale][flask][isolate]):
                    num_cols += 1
                    print(pair, ale, flask, isolate, replicate)
                    read_coverage = coverage_dict[ale][flask][isolate][replicate]
                    summary_file = return_summary_file_loc(ale, flask,
                                                           isolate, replicate)
                    mean = get_fit_mean_from_breseq(summary_file)

                    read_cov = []
                    for i, cov in enumerate(read_coverage):
                        filtered_cov = cov if starting_dup_pos[i] == 0 else mean
                        read_cov.append(filtered_cov)

                    sectioned_cov, _ = \
                        get_average_and_stdev_of_sections(read_cov, genome_size,
                                                          num_sections=sections)
                    y = []
                    for value in sectioned_cov:
                        if value > (4 * mean) and not unfiltered:
                            value = 4 * mean
                        y.append(value)

                    x = np.array(range(0, len(y))) * genome_size / len(y)

                    xy = lowess(y, x, it=0, frac=.005)
                    ymax = max(xy[:, 1] / mean) if max(
                        xy[:, 1] / mean) > 2 else 2

                    # a rows num_cols columns
                    if isolate not in [30, 40]:
                        ax = axes[num_cols][a]
                        kind = ''
                    elif isolate == 30:
                        ax = axes[-2][a]
                        kind = 'endclone'
                    elif isolate == 40:
                        ax = axes[-1][a]
                        kind = 'endclone'
                    add_gene_positions(ax, ymax, kind='all')

                    plot_and_format(ax, xy[:, 0], xy[:, 1], genome_size, ale,
                                    flask, isolate, replicate, mean)

                    ax.plot([0, genome_size], [cutoff, cutoff], '--r',
                            label='cutoff', linewidth=3)

                    if ale == 5:
                        ax_gltJ = axes_gltJ[0][num_cols]
                        plot_and_format(ax_gltJ, xy[:, 0], xy[:, 1],
                                        genome_size, ale, flask, isolate,
                                        replicate, mean, xrange=[5e5, 15e5])
                        add_gene_positions(ax_gltJ, ymax, kind='gltJ')
                        ax_gltJ.set_title(ax_gltJ.get_title().replace('Ale', '$\Delta$hisD & $\Delta$gdhAgltB\nAle'))
                        ax_gltJ.set_ylim([0, 25])
                        ax_gltJ.ticklabel_format(style='sci', axis='x',
                                                 scilimits=(0, 0))

                    elif ale == 11:
                        ax_gltJ = axes_gltJ[1][num_cols]
                        plot_and_format(ax_gltJ, xy[:, 0], xy[:, 1],
                                        genome_size, ale, flask, isolate,
                                        replicate, mean, xrange=[2e6, 3.3e6])
                        add_gene_positions(ax_gltJ, ymax, kind='hisJ')
                        ax_gltJ.set_title(ax_gltJ.get_title().replace('A', '$\Delta$hisD & $\Delta$gltAprpC\nA'))
                        ax_gltJ.set_ylim([0, 5])
                        ax_gltJ.ticklabel_format(style='sci', axis='x',
                                                 scilimits=(0, 0))

    save_fig(fig, pair, '', unfiltered=unfiltered)
    save_fig(fig_gltJ, pair, 'gltJ', unfiltered=unfiltered)
