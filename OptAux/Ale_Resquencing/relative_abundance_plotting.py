from __future__ import absolute_import, division, print_function

from os.path import dirname, abspath
from collections import Counter

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

pwd = dirname(abspath(__file__))
colors = sns.color_palette(palette='muted')

plt.rcParams['font.size'] = 15.
plt.rcParams['font.weight'] = 'bold'

plt.rcParams['legend.fontsize'] = 15.0
plt.rcParams['legend.framealpha'] = .6
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['grid.color'] = 'white'
plt.rcParams['figure.edgecolor'] = 'black'
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['axes.grid'] = False

strain_to_ko = {'hisD': r'$\Delta \mathrm{hisD}$',
                'pyrC': r'$\Delta \mathrm{pyrC}$',
                'gltA': r'$\Delta \mathrm{gltA} \Delta \mathrm{prpC}$',
                'gltB': r'$\Delta \mathrm{gdhA} \Delta \mathrm{gltB}$'}


def _get_x_label_array(df):
    out = [''] * len(df.columns)
    ales = []
    for column in df.columns:
        ales.append(column[0])
    counted_ales = Counter(ales)
    i = 0
    for ale in sorted(counted_ales):
        if ale % 1 != 0 or ale == 7:
            i += 1
        else:
            out[(i + counted_ales[ale]//2)] = '%i' % ale
            i += int(counted_ales[ale])
    return out


def _add_space_between_ales(df):

    df = df.T.reset_index()
    last_ale = df['Ale'][0]
    i = 0
    start_number = len(df.index)

    for ale in list(df['Ale']):
        if ale != last_ale:
            new_ale = (ale + last_ale) / 2
            df.loc[start_number + i, 'Ale'] = new_ale
            df.loc[start_number + i, 'Flask'] = 0
            df.loc[start_number + i, 'Isolate'] = 0
            i += 1
        last_ale = ale

    df = df.set_index(['Ale', 'Flask', 'Isolate'])
    df = df.sort_index(axis=0)

    return df.T


def plot_abundance_df(abundance_df, strain_1, strain_2, type, ax):

    abundance_df = _add_space_between_ales(abundance_df)
    left = np.array(list(range(len(abundance_df.columns))))

    ax.bar(left, abundance_df.loc[strain_1].values, label=strain_to_ko[strain_1])
    ax.bar(left, abundance_df.loc[strain_2].values,
           bottom=abundance_df.loc[strain_1].values, label=strain_to_ko[strain_2])

    ax.set_xticks(left)
    label_array = _get_x_label_array(abundance_df)
    ax.set_xticklabels(label_array, ha='center')
    ax.plot([-.5, len(left) - .5], [1, 1], 'r--')
    ax.set_xlim([-.5, len(left)-.5])
    ax.set_title(('Abundance by %s' % type).replace('characteristic', '\n characteristic mutations'))
    ax.legend(loc='lower left', facecolor='w', framealpha=.8, frameon=True)


def plot_pairwise_comparison_of_preditions(table_loc, save_loc,
                                           normalize=True):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    for i, pair in enumerate(['hisD_gltA', 'hisD_gltB', 'hisD_pyrC']):
        s1, s2 = pair.split('_')
        ax = axes[i]

        # Upload dataframes of abundance predictions based on coverage and
        # characteristic mutations (unnormalized)
        df_cov = pd.read_csv(
            '%s/abundance_by_coverage_%s.csv' % (
                table_loc, pair.replace('_', '')),
            index_col=0)
        df_char = pd.read_csv(
            '%s/abundance_by_characteristic_%s.csv' % (
                table_loc, pair.replace('_', '')),
            index_col=0)

        if normalize:
            def _normalize_values(df):
                s1_unnormalized = np.array(df.T[s1].values)
                s2_unnormalized = np.array(df.T[s2].values)
                df.T[s1] = df.T[s1] / (s1_unnormalized + s2_unnormalized)
                df.T[s2] = df.T[s2] / (s1_unnormalized + s2_unnormalized)
                return df
            df_cov = _normalize_values(df_cov)
            df_char = _normalize_values(df_char)

        # Append mutation identifiers with prediction type and merge dfs
        df_cov = df_cov.rename_axis(lambda x: x + '_coverage')
        df_char = df_char.rename_axis(lambda x: x + '_characteristic')
        df = df_char.append(df_cov).T

        ax.scatter(df_char.T['%s_characteristic' % s1],
                   df_cov.T['%s_coverage' % s1], label=strain_to_ko[s1])
        ax.scatter(df_char.T['%s_characteristic' % s2],
                   df_cov.T['%s_coverage' % s2], label=strain_to_ko[s2])

        # Calculate Correlation
        cor1 = df[['%s_characteristic' % s1, '%s_coverage' % s1]].corr().iloc[
            0, 1]
        cor2 = df[['%s_characteristic' % s2, '%s_coverage' % s2]].corr().iloc[
            0, 1]
        ax.plot([0, 1], [0, 1])
        ax.legend()
        ax.set_xlabel('Characteristic Mutations')

    axes[1].set_title('Comparison of Methods')
    axes[0].set_ylabel('Gene Deletion Coverage')
    fig.tight_layout()
    fig.savefig('%s/cov_char_comparison.png' % save_loc)
