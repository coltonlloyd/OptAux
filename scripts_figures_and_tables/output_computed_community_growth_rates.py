from __future__ import print_function, absolute_import, division

from glob import glob
import os
import json
import tarfile
from collections import defaultdict
from textwrap import wrap

import pandas
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from cycler import cycler

import cobra
from optaux import resources
from optaux.resources.update_mpl_rcparams import update_rcparams

update_rcparams()
colors = ['#0099E6', '#F23814', '#DFB78B', '#9E9E9E', '#12255A']
colors.extend(reversed(colors[:-1]))
plt.rcParams['axes.prop_cycle'] = cycler(color=colors)

here = os.path.abspath(os.path.dirname(__file__))
resource_dir = resources.__path__[0]
m_model_loc = resource_dir + '/iJO1366.json'
ijo = cobra.io.load_json_model(m_model_loc)

cov_abundance_df = \
    pandas.read_excel('%s/relative_abundance/tables/coverage_abundances_strain1.xlsx' %
                      here, index_col=0, header=None)

rxn_to_gene = {'CS': r'$\Delta \mathit{gltA} \Delta \mathit{prpC}$',
               'HISTD': r'$\Delta \mathit{hisD}$',
               'GLUDy:GLUSy': r'$\Delta \mathit{gdhA} \Delta \mathit{gltB}$',
               'DHORTS': r'$\Delta \mathit{pyrC}$'}

legend_order = {'HISTD-CS': ['glu__L', 'gln__L', 'pro__L', 'orn'],
                'HISTD-DHORTS': ['orot', 'csn'],
                'HISTD-GLUDy:GLUSy': ['gln__L', 'ala__L', 'asn__L', 'orn',
                                      'arg__L']}


def remove_strain_from_id(s):
    return s.replace('_S1', '').replace('_S2', '')


def add_experimental_abundance_boxplot(abundance_df, ale, axis):
    ale_name = '__'.join([i.split(':')[0] for i in ale.split('-')])
    if ale_name in abundance_df.index:
        abundance = abundance_df.loc[ale_name].dropna().values
        sns.boxplot(data=abundance, orient='horizontal', ax=axis)
        sns.swarmplot(data=abundance, size=5, color="0", linewidth=0, ax=axis,
                      orient='horizontal')
        axis.set_yticks([0])
        axis.set_yticklabels(['Experimentally \n Inferred'], rotation=90,
                             va='center', ha='center')
        axis.tick_params(axis='y', pad=15)
        axis.set_xlabel(r'Fraction $\Delta \mathrm{hisD}$')


def filter_based_on_plot_kind(plot_kind, unmodeled_name, pair, line_num=0):
    line_kwargs = {}
    if plot_kind == 'metabolite_limitation':
        suffix = unmodeled_name.replace('EX_', '').replace('_e', '')
        save_prefix = 'met'

        line_kwargs['marker'] = 'o'
        line_kwargs['linewidth'] = 3.
        line_kwargs['markersize'] = 10
    elif plot_kind == 'unmodeled_sweep':
        q1, q2 = unmodeled_name.split('_')[:2]
        suffix = r'%s - %s' % (q1 + '%', q2 + '%')
        save_prefix = 'unmodeled'

        if float(q2) > float(q1):
            line_kwargs['linestyle'] = ':'
            line_kwargs['linewidth'] = 2.
            line_kwargs['marker'] = 'o'
        elif float(q2) < float(q1):
            line_kwargs['linestyle'] = '--'
            line_kwargs['linewidth'] = 2.
            line_kwargs['marker'] = '^'
        elif float(q2) == float(q1):
            line_kwargs['linestyle'] = '-'
            line_kwargs['linewidth'] = 6.
    elif plot_kind == 'secretion_keff_sweep':
        k1, k2 = unmodeled_name.split('_')[:2]
        k1 = float(k1)
        k2 = float(k2)
        if (k1 / k2) < 1.:
            suffix = '%.4f' % (k1 / k2)
            suffix = suffix.rstrip('0')
        else:
            suffix = '%i' % (k1 / k2)
        save_prefix = 'secretion_keff'
        if float(k2) > float(k1):
            line_kwargs['linestyle'] = '-'
            line_kwargs['linewidth'] = 4.
        elif float(k2) < float(k1):
            line_kwargs['linestyle'] = ':'
            line_kwargs['linewidth'] = 4.
        elif float(k2) == float(k1):
            line_kwargs['marker'] = 'o'
            line_kwargs['markersize'] = 10

    elif plot_kind == 'default':
        suffix = ''
        save_prefix = 'default'
        line_kwargs['linestyle'] = '-'
        line_kwargs['marker'] = 'o'
    else:
        raise UserWarning('Not valid plot_kind')

    return suffix, save_prefix, line_kwargs


def get_growth_rates(sim_loc, ale, unmodeled, me_model=True):

    sims = glob('/'.join([sim_loc, ale, unmodeled, '*_sol.json']))
    if len(sims) == 0:
        return None, None
    me_gr = []
    x = []
    for sim in sorted(sims):
        # fraction of strain 1
        frac = (float(sim.split('/')[-1].split('_')[0]))
        with open(sim, 'r') as f:
            x_dict = json.load(f)
        x.append(frac)
        if me_model:
            me_gr.append(x_dict['biomass_dilution_S1'])
        else:
            me_gr.append(x_dict['biomass'])

    return x, me_gr


def make_default_comparison_plot(abundance_df, me_sim_dir, m_sim_dir):
    keff_sets = ['65', 'default', 'null', 'ML']
    me_plot_kinds = ['default']
    gridkw = dict(height_ratios=[5, 1, 1])
    fig1, (axes1, axes2, axes3) = plt.subplots(3, 3, sharex='col',
                                        gridspec_kw=gridkw,
                                        figsize=(15, 6))

    ale_plot_locs = {}
    max_gr_positions = defaultdict(list)
    for keff_set in keff_sets:
        for plot_kind in me_plot_kinds:
            # Glucose limited plots are identical
            if keff_set == 'default' and plot_kind == 'glucose_limited':
                continue
            sim_loc = '/'.join(['/home/sbrg-cjlloyd/Desktop/community_sims_output_%s_keffs' % keff_set, plot_kind])
            for i, directory in enumerate(glob('/'.join([sim_loc, '*']))):
                ale = directory.split('/')[-1]
                if ale not in ale_plot_locs:
                    ale_plot_locs[ale] = i
                else:
                    i = ale_plot_locs[ale]
                annotation = [x.split('/')[-1] for
                              x in glob('/'.join([sim_loc, ale, '*']))][0]

                x, gr = get_growth_rates(sim_loc, ale, annotation)
                if not x:
                    continue
                x = [0] + x + [1]
                gr = [0] + gr + [0]
                gr_array = np.array(gr)
                norm_gr = gr_array/gr_array.max()
                axes1[i].plot(x, norm_gr, label=plot_kind+'_'+keff_set,
                              linewidth=3)

                x_max = x[np.argmax(gr_array)]
                x_array = np.array(x)
                x_close_to_max =x_array[norm_gr > .95]
                gr_close_to_max = norm_gr[norm_gr > .95]
                max_gr_positions[ale].extend(x_close_to_max)
                gr_max = gr_array.max()
                axes1[i].arrow(x_max, .1, 0, -.1, width=.01,
                               color=axes1[i].get_lines()[-1].get_color())
                axes1[i].fill_between(x_close_to_max, gr_close_to_max,
                                      np.zeros(len(x_close_to_max)),
                                      facecolor='#B6C3C5',
                                      alpha=.2
                                      )
    for ale, position in ale_plot_locs.items():
        add_experimental_abundance_boxplot(abundance_df, ale, axes3[position])

        x, gr = get_growth_rates(m_sim_dir, ale, 'no_unmodeled', me_model=False)
        x = [0] + x + [1]
        gr = [0] + gr + [0]
        norm_gr = np.array(gr) / np.array(gr).max()
        #axes1[position].plot(x, norm_gr, label='m_model',
        #                      linewidth=3)
        print(max_gr_positions[ale])
        axes2[position].set_xlim([0, 1])
        #sns.boxplot(np.array(max_gr_positions[ale]), ax=axes2[position])
        #sns.swarmplot(max_gr_positions[ale], size=5, color="0", linewidth=0,
        #             ax=axes2[position], orient='horizontal')
        bins = (np.array(max_gr_positions[ale]).max() - \
               np.array(max_gr_positions[ale]).min()) * 10
        print(np.array(max_gr_positions[ale]).max())
        print(np.array(max_gr_positions[ale]).min())
        print(int(bins))
        sns.distplot(np.array(max_gr_positions[ale]), ax=axes2[position], rug=True,
                     bins=int(bins)+1)
        axes1[position].legend()
        axes1[position].set_title(ale)
    return fig1


def make_computational_abundance_plots(abundance_df, sim_dir, plot_kind,
                                       normalize_to_max, max_gr_dict,
                                       subplot_axes=None):
    sim_loc = '/'.join([sim_dir, plot_kind])
    if plot_kind == 'default' or plot_kind == 'unmodeled_sweep':
        plot_experimental_abundance = True
        subplot_axes = None
        plot = True
    else:
        plot = False
        plot_experimental_abundance = False

    if plot_experimental_abundance and subplot_axes is None:
        gridkw = dict(height_ratios=[5, 1])
        fig1, (axes1, axes2) = plt.subplots(2, 3, sharex='col',
                                            gridspec_kw=gridkw,
                                            figsize=(15, 5))
    elif subplot_axes is None:
        fig1, axes1 = plt.subplots(1, 3, sharex='row', sharey='col',
                                   figsize=(15, 4))
    else:
        axes1 = subplot_axes
        axes1[0].set_ylabel('Computed Community \n'
                            r'Growth Rate ($hr^{-1}$)')
    # Plot growth rates, reversed to get HISTD-GLUDy:GLUSy first
    for i, ale in enumerate(glob('/'.join([sim_loc, '*']))):
        ale = ale.split('/')[-1]

        print(ale)
        if ale not in max_gr_dict:
            max_gr_dict[ale] = defaultdict(list)

        ax1 = axes1[i]
        if plot_experimental_abundance:
            ax2 = axes2[i]
            add_experimental_abundance_boxplot(abundance_df, ale, ax2)

        unmodeled_list = [i.split('/')[-1] for
                          i in glob('/'.join([sim_loc, ale, '*']))]

        # [(strain_1 unmodeled/keff multiplier, filename)]
        if plot_kind != 'metabolite_limitation':
            number_list = [(float(i.split('_')[0]), i) for i in unmodeled_list]
        else:
            number_list = [(i, k) for i, k in enumerate(sorted(unmodeled_list))]

        # Get dictionary of optimal hisD fraction to metabolites w/ that
        # optimal fractoin
        max_gr_frac_to_met = defaultdict(list)
        for i, unmodeled in sorted(number_list):
            x, me_gr = get_growth_rates(sim_loc, ale, unmodeled)
            if not x:
                continue
            max_gr_frac_to_met[x[np.argmax(me_gr)]].append(
                unmodeled.replace('EX_', '').replace('_e', ''))
        max_gr_frac_to_color = {k: i for i, k in enumerate(max_gr_frac_to_met)}

        for i, unmodeled in sorted(number_list):
            print(unmodeled)

            out = filter_based_on_plot_kind(plot_kind, unmodeled, ale,
                                            line_num=int(i))
            if out is True:
                continue
            else:
                suffix, save_prefix, line_kwargs = out

            x, me_gr = get_growth_rates(sim_loc, ale, unmodeled)
            # if None is returned for x, continue
            if not x:
                continue

            # store values in a dictionary for box plots
            me_max_gr = np.array(me_gr).max()
            optimal_fraction = x[np.argmax(me_gr)]
            if plot_kind != 'metabolite_limitation':
                mult1, mult2 = unmodeled.split('_')[:2]

                # take inverse to make bar plot in same order
                if plot_kind == 'secretion_keff_sweep':
                    mult1 = 1 / float(mult1)
                    mult2 = 1 / float(mult2)
                else:
                    mult1 = float(mult1)
                    mult2 = float(mult2)
                if mult1 > mult2:
                    max_gr_dict[ale][plot_kind + '_1'].append(optimal_fraction)
                elif mult1 < mult2:
                    max_gr_dict[ale][plot_kind + '_2'].append(optimal_fraction)
            else:
                max_gr_dict[ale][plot_kind].append(optimal_fraction)
                print('max', unmodeled, me_max_gr)

            if normalize_to_max:
                me_gr = np.array(me_gr) / me_max_gr

            # growth rate of zero at 0% and 100% strain 1 fraction
            if plot_kind == 'metabolite_limitation':
                if optimal_fraction in max_gr_frac_to_met:
                    legend = ', '.join(max_gr_frac_to_met.pop(optimal_fraction))
                    wrap_legend = '\n'.join(wrap(legend, 40))
                    ax1.plot([0] + x + [1], [0] + me_gr + [0],
                             label='%s' % (wrap_legend),
                             color=colors[max_gr_frac_to_color[optimal_fraction]],
                             **line_kwargs)
                else:
                    ax1.plot([0] + x + [1], [0] + me_gr + [0],
                             color=colors[
                                       max_gr_frac_to_color[optimal_fraction]],
                                   **line_kwargs)
            else:
                ax1.plot([0] + x + [1], [0] + me_gr + [0],
                         label='%s' % suffix, **line_kwargs)

            # grey region of max growth on plot
            if plot_kind == 'default':
                me_max_loc = np.argmax(me_gr)
                ymin, ymax = ax1.get_ylim()
                ax1.set_facecolor('white')
                ax1.fill_between([x[me_max_loc-1], x[me_max_loc],
                                  x[me_max_loc+1]],
                                 [me_gr[me_max_loc-1], me_gr[me_max_loc],
                                  me_gr[me_max_loc+1]],
                                 [ymin, ymin, ymin], facecolor='#B6C3C5',
                                 alpha=.3)
                ax1.plot([x[me_max_loc], x[me_max_loc]],
                         [ymin, me_gr[me_max_loc]], 'k--')
                # plotting changes ymin, set again
                ax1.set_ylim([ymin, ymax])

            # Add legend for both plot types
            if plot_kind != 'metabolite_limitation':
                lgd = axes1[1].legend(ncol=2, frameon=True,
                                      title=r'$\mathbf{\frac{K_{eff} \ Export '
                                            r'\ \Delta hisD}{K_{eff} \ Export '
                                            r'\ Partner \ Strain}}$')
                # Only middle plot has legend
                if lgd:
                    plt.setp(lgd.get_title(), fontsize=15)
            else:
                ax1.legend()
                handles, labels = ax1.get_legend_handles_labels()
                order_list = legend_order[ale]
                sorted_handles = [''] * len(order_list)
                sorted_labels = [''] * len(order_list)

                for i, met in enumerate(order_list):
                    for h, l in zip(handles, labels):
                        if met in l:
                            sorted_handles[i] = h
                            sorted_labels[i] = l
                lgd = ax1.legend(handles, labels, ncol=1,
                                 frameon=True,
                                 title=r'$\mathbf{Metabolite \ '
                                       r'Cross\mathrm{-}fed}$')
                if lgd:
                    plt.setp(lgd.get_title(), fontsize=12)

        if plot:
            axes1[0].set_ylabel('Computed Community \n'
                                r'Growth Rate ($hr^{-1}$)')
            ax1.set_xlim([0, 1.])
            ax1.set_title(' & '.join([rxn_to_gene[i] for i in ale.split('-')]))
            fig1.tight_layout()
            fig1.savefig('%s/community_plots/%s/%s.png' % (
                         here, plot_kind, save_prefix))


def make_sweeps_box_plots(max_gr_dict, abundance_df, plot_kind,
                          subplot_axes=None):
    if subplot_axes is None:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey='row')
    else:
        axes = subplot_axes

    for i, ale in enumerate(['HISTD-CS', 'HISTD-DHORTS', 'HISTD-GLUDy:GLUSy']):
        ax = axes[i]
        ale_name = '__'.join([i.split(':')[0] for i in ale.split('-')])
        sns.boxplot(data=[max_gr_dict[ale]['metabolite_limitation'],
                          max_gr_dict[ale]['%s_1' % plot_kind],
                          max_gr_dict[ale]['%s_2' % plot_kind],
                          abundance_df.loc[ale_name].dropna().values,
                          cov_abundance_df.loc[ale_name].dropna().values],
                    orient='horizontal', ax=ax)
        sns.swarmplot(data=[max_gr_dict[ale]['metabolite_limitation'],
                            max_gr_dict[ale]['%s_1' % plot_kind],
                            max_gr_dict[ale]['%s_2' % plot_kind],
                            abundance_df.loc[ale_name].dropna().values,
                            cov_abundance_df.loc[ale_name].dropna().values],
                      size=5, color="0", linewidth=0, ax=ax,
                      orient='horizontal')
        ax.set_yticks([0, 1, 2, 3, 4])

        ax.set_yticklabels(["Metabolite restriction \n" +
                            r"($\mathbf{Panel \ A}$)",

                            r"$\Delta$hisD" + " export less efficient \n " +
                            r"($\mathbf{Panel \ B, \ solid}$)",

                            r"$\Delta$hisD" + " export more efficient \n " +
                            r"($\mathbf{Panel \ B, \ dashed}$)",

                            "Experimentally Inferred \n " +
                            r"($\mathbf{By \ mutation}$)",

                            "Experimentally Inferred \n " +
                            r"($\mathbf{By \ coverage}$)"],

                           va='center', ha='center')
        ax.tick_params(axis='y', which='major', pad=90)
        ax.set_xlim([0, 1.])
        # ax.set_title(' & '.join([rxn_to_gene[i] for i in ale.split('-')]))
    axes[1].set_xlabel(r"$\Delta \mathit{hisD}$ Fraction", fontsize=23,
                       color='black')


if __name__ == '__main__':
    abundance_df = \
        pandas.read_excel('%s/relative_abundance/tables/'
                          'characteristic_abundances_strain1.xlsx' % here,
                          index_col=0, header=None)

    plot_kinds = ['metabolite_limitation', 'secretion_keff_sweep']
    normalize_to_max = False
    main_text_figure = True

    sim_dir = '%s/community_sims_output_ML_keffs/' % '/home/sbrg-cjlloyd/Desktop'
    save_location = '%s/community_plots/' % here
    m_sim_dir = '%s/community_m_sims/' % '/home/sbrg-cjlloyd/Desktop'

    # ######################### Make M + ME comparison
    plt.rcParams['axes.facecolor'] = 'w'
    fig = make_default_comparison_plot(abundance_df, sim_dir, m_sim_dir)
    fig.savefig('/home/sbrg-cjlloyd/Desktop/test.png')
    # ######################### Make Figure 8 ################################
    max_gr_dict = {}
    fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharey='row',
                             sharex='col')
    for i, plot_kind in enumerate(plot_kinds):

        # Get location of simulation fluxes based on simulation type
        make_computational_abundance_plots(abundance_df, sim_dir, plot_kind,
                                           normalize_to_max, max_gr_dict,
                                           subplot_axes=axes[i])

    make_sweeps_box_plots(max_gr_dict, abundance_df, plot_kind,
                          subplot_axes=axes[2])

    # Add titles to top row
    for i, ale in enumerate(['HISTD-CS', 'HISTD-DHORTS', 'HISTD-GLUDy:GLUSy']):
        axes[0][i].set_title(' & '.join([rxn_to_gene[i] for
                                         i in ale.split('-')]),
                             fontsize=25)

    # Add labels to each panel
    for i, t in enumerate(['A', 'B', 'C']):
        fig.text(0.085, 1-i/3.2, t, fontsize=35, fontweight='bold',
                 va='top', color='black')

    fig.tight_layout()
    fig.savefig('%s/community_plots/Figure_8.png' % here,  dpi=250)

    # ################### Make supplementary fig ##############################
    plt.rcParams['axes.facecolor'] = 'w'
    max_gr_dict = {}
    plot_kind = 'default'

    save_loc = save_location + plot_kind
    if not os.path.exists(save_loc):
        os.mkdir(save_loc)

    fig, ax = plt.subplots(1, 1)
    make_computational_abundance_plots(abundance_df, sim_dir, plot_kind,
                                       normalize_to_max, max_gr_dict,
                                       subplot_axes=ax)
    fig.savefig('%s/community_plots/Supplementary_fig.png' % here, dpi=250)
