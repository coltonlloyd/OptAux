from __future__ import print_function, absolute_import, division
import matplotlib
matplotlib.use('agg')

from glob import glob
import pickle
import pandas
import numpy as np
import os
import json

from matplotlib import pyplot as plt

from os.path import expanduser
from collections import defaultdict
import seaborn as sns
from cycler import cycler
from colour import Color
from IPython import embed
from textwrap import wrap

import cobra
from optaux.me_community.me_model_community import (make_binary_community_me,
                                                    scale_exchange_fluxes)
from optaux.me_community.run_ALE_pairs import restrict_uptake
from optaux import resources
from optaux.resources import possible_uptake

here = os.path.abspath(os.path.dirname(__file__))
resource_dir = resources.__path__[0]
m_model_loc = resource_dir + '/iJO1366.json'
ijo = cobra.io.load_json_model(m_model_loc)

skip_list = ['h2o', 'etoh', 'o2', 'mg2', 'nh4', 'so4', 'h', 'ac', 'for',
             'pi', 'lac__D', 'co2', 'glc__D', 'fe2', 'fe3', 'fum']

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.facecolor'] = '#eeeeee'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['text.color'] = 'black'
plt.rcParams['grid.color'] = 'white'
plt.rcParams['figure.edgecolor'] = 'black'
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['axes.grid'] = False
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['legend.edgecolor'] = 'k'
plt.rcParams['legend.framealpha'] = .8
plt.rcParams['legend.fancybox'] = True

red = Color('#0099E6')

colors = ['#0099E6', '#F23814', '#DFB78B', '#9E9E9E', '#12255A']
colors.extend(reversed(colors[:-1]))
plt.rcParams['axes.prop_cycle'] = cycler(color=colors)

#[i.get_hex() for i in list(red.range_to(Color('#12255A'), 10))]

red_colors = list(red.range_to(Color('#12255A'), 10))

cov_abundance_df = \
    pandas.read_excel('%s/relative_abundance/tables/coverage_abundances_strain1.xlsx' %
                      here, index_col=0, header=None)

rxn_to_gene = {'CS': r'$\Delta gltAprpC$',
               'HISTD': r'$\Delta hisD$',
               'GLUDy:GLUSy': r'$\Delta gdhAgltB$',
               'DHORTS': r'$\Delta pyrC$'}

legend_order = {'HISTD-CS': ['glu__L', 'gln__L', 'pro__L', 'orn'],
                'HISTD-DHORTS': ['orot', 'csn'],
                'HISTD-GLUDy:GLUSy': ['gln__L', 'ala__L', 'asn__L', 'orn', 'arg__L']}


def remove_strain_from_id(s):
    return s.replace('_S1', '').replace('_S2', '')


def get_metabolite_exchange(x_dict, exchange_fluxes, frac_S1):

    sim_fluxes = {}
    for r, v in x_dict.items():

        if 'reverse' in r or not r.startswith('EX_') or '_Shared' in r:
            continue

        met_id = remove_strain_from_id(r).replace('EX_', '').replace('_e', '')
        if met_id in skip_list:
            continue

        sim_fluxes[r] = v * frac_S1 if '_S1' in r else -v * (1-frac_S1)

    skip_reaction = []
    for r in list(sim_fluxes.keys()):
        if r in skip_reaction:
            continue
        skip_reaction.append(r)
        v = sim_fluxes[r]
        s2_r = r.replace('_S1', '_S2')
        s1_r = r.replace('_S2', '_S1')
        if '_S1' in r and s2_r not in skip_reaction:
            v += sim_fluxes.get(s2_r, 0)
            skip_reaction.append(s2_r)
        elif '_S2' in r and s1_r not in skip_reaction:
            v += sim_fluxes.get(s1_r, 0)
            skip_reaction.append(s1_r)

        exchange_fluxes[remove_strain_from_id(r)].append(v)

    return exchange_fluxes


def filter_exchange_fluxes(exchange_fluxes, threshold=.2):
    output_fluxes = {}
    for r, flux in exchange_fluxes.items():

        flux = np.array(flux)

        if len(flux) == 0 or abs(flux).sum() <= threshold:
            continue

        output_fluxes[r] = flux

    return output_fluxes


def make_community_m_model(ijo, fraction_strain1, ko1s, ko2s,
                           scale_uptake=False):
    ijo1 = ijo.copy()
    ijo2 = ijo.copy()

    com = make_binary_community_me(ijo1, ijo2, ME_model=False)

    # add combined biomass reaction
    r = com.reactions.BIOMASS_Ec_iJO1366_core_53p95M_S1 + \
        com.reactions.BIOMASS_Ec_iJO1366_core_53p95M_S2
    r.id = 'COMMUNITY_BIOMASS'
    com.add_reaction(r)
    com.objective = r.id

    com.reactions.EX_glc__D_e_Shared.lower_bound = -5.6

    for ko in ko1s:
        com.reactions.get_by_id(ko + '_S1').knock_out()
    for ko in ko2s:
        com.reactions.get_by_id(ko + '_S2').knock_out()

    # Make sure strains cannot secrete any metabolites that it is auxotorphic
    # for
    for r in possible_uptake.ko_uptakes[str(ko1s)]:
        com.reactions.get_by_id(r + '_S1').upper_bound = 0
        print(str(ko1s), r)
    for r in possible_uptake.ko_uptakes[str(ko2s)]:
        com.reactions.get_by_id(r + '_S2').upper_bound = 0
        print(str(ko2s), r)

    restrict_uptake(com, ko1s, ko2s, is_me_model=True)

    scale_exchange_fluxes(com, fraction_strain1, uptake=scale_uptake)

    return com


def make_exchange_met_plots(ex, x, unmodeled, ale, m_model=True, axis=None,
                            num=0):
    if not axis:
        fig2, ax = plt.subplots()
        plot=True
    else:
        ax = axis
        plot=False

    i = 0
    for r, flux in ex.items():
        if i < 5:
            line_type = '-'
        elif i < 10:
            line_type = '--'
        elif i < 15:
            line_type = '-.'
        else:
            line_type = ':'

        met = r.replace('EX_', '').replace('_e', '')
        color = colors[0] if met == 'his__L' else colors[1]
        ax.plot(x, flux, line_type, linewidth=5,
                label=met, color=color)
        i += 1

    ax.legend(ncol=1, fontsize=15)
    ax.set_title(ale)
    ax.set_xlim([.05, .95])
    if plot:
        ax.legend(ncol=2, fontsize=15)
        ax.set_xlabel('Fraction Strain 1')
        ax.set_ylabel('Metabolite Exchange Rate')
        if m_model:
            fig2.savefig(expanduser('%s/community_plots/m_model_%s_%s.png') % (
                here, unmodeled, ale))
        else:
            fig2.savefig(expanduser('%s/community_plots/me_model_%s_%s.png') % (
                here, unmodeled, ale))


def add_experimental_abundance_boxplot(abundance_df, ale, axis):
    ale_name = '__'.join([i.split(':')[0] for i in ale.split('-')])
    if ale_name in abundance_df.index:
        abundance = abundance_df.loc[ale_name].dropna().values
        sns.boxplot(data=abundance, orient='horizontal', ax=axis)
        axis.set_yticks([0])
        axis.set_yticklabels(['Experiment'], rotation=90, va='center')
        axis.set_xlabel('Fraction Strain 1')


def filter_based_on_plot_kind(plot_kind, unmodeled_name, pair, line_num=0):
    line_kwargs = {}
    if plot_kind == 'metabolite_limitation':
        suffix = unmodeled_name.replace('EX_', '').replace('_e', '')
        save_prefix = 'met'
        if line_num < 5:
            line_kwargs['marker'] = 'o'
        elif line_num < 10:
            line_kwargs['marker'] = 'o'
        elif line_num < 15:
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


def get_growth_rates(sim_loc, ale, unmodeled):
    output_dir = '/'.join([here, 'community_me_growth_rates', ale, unmodeled])
    if os.path.exists('%s/all_grs.json' % output_dir):
        with open('%s/all_grs.json' % output_dir, 'r') as f:
            x, me_gr = json.load(f)
    else:
        sims = glob('/'.join([sim_loc, ale, unmodeled, '*.pickle']))
        if len(sims) == 0:
            return None, None
        me_gr = []
        x = []
        for sim in sorted(sims):
            # fraction of strain 1
            frac = (float(sim.split('/')[-1].split('_')[0]))
            if int(100. * frac) % 5 != 0 or int(100. * frac) % 10 == 0:
                continue
            print(sim)
            x.append(frac)
            with open(sim, 'rb') as f:
                me_sol = pickle.load(f, encoding='latin1')

            me_gr.append(me_sol.x_dict['biomass_dilution_S1'])

        os.makedirs(output_dir)
        with open('%s/all_grs.json' % output_dir, 'w') as f:
            json.dump([x, me_gr], f)

    return x, me_gr


def get_me_exchange_array(sim_loc, ale, unmodeled):
    me_ex = defaultdict(list)
    x = []
    sims = glob('/'.join([sim_loc, ale, unmodeled, '*.pickle']))
    if len(sims) == 0:
        return None, None

    for sim in sorted(sims):
        # fraction of strain 1
        frac = (float(sim.split('/')[-1].split('_')[0]))
        if int(100. * frac) % 5 != 0 or int(100. * frac) % 10 == 0:
            continue
        print(sim)
        x.append(frac)
        with open(sim, 'rb') as f:
            me_sol = pickle.load(f, encoding='latin1')

        me_ex = get_metabolite_exchange(me_sol.x_dict, me_ex, frac)

    return x, me_ex


def make_computational_abundance_plots(abundance_df, sim_loc, plot_kind,
                                       normalize_to_max, max_gr_dict,
                                       plot_exchange=False, subplot_axes=None):

    if plot_kind == 'default':
        plot_experimental_abundance = True
        plot_exchange = True
        subplot_axes = None
        plot = True
    else:
        plot_experimental_abundance = False

    plot_exchange = True

    if plot_kind in ['default', 'secretion_keff_sweep']:
        exchange_fig, exchange_axes = \
                plt.subplots(1, 3, sharex=True, figsize=(15, 4))
        single_exchange_plot = True
    else:
        single_exchange_plot = False
        exchange_axes = None

    if plot_experimental_abundance and subplot_axes is None:
        gridkw = dict(height_ratios=[5, 1])
        fig1, (axes1, axes2) = plt.subplots(2, 3, sharex=True, gridspec_kw=gridkw,
                                            figsize=(15, 5))
    elif subplot_axes is None:
        fig1, axes1 = plt.subplots(1, 3, sharex=True, sharey=True,
                                            figsize=(15, 4))
    else:
        axes1 = subplot_axes
        axes1[0].set_ylabel(r'Community Growth Rate ($hr^{-1}$)')
    # Plot growth rates, reversed to get HISTD-GLUDy:GLUSy first
    for i, ale in enumerate(glob('/'.join([sim_loc, '*']))):
        if "HISTD" not in ale:
            continue
        ale = ale.split('/')[-1]
        plot = True
        print(ale)
        if ale not in max_gr_dict:
            max_gr_dict[ale] = defaultdict(list)

        ax1 = axes1[i]
        if single_exchange_plot:
            exchange_ax = exchange_axes[i]

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
                                            line_num=i)
            if out is True:
                continue
            else:
                suffix, save_prefix, line_kwargs = out

            x, me_gr = get_growth_rates(sim_loc, ale, unmodeled)
            # if None is returned for x, continue
            if not x:
                continue

            if plot_exchange:
                print(unmodeled)
                if single_exchange_plot and unmodeled == '100.00_100.00_0_secretion_multiplier':
                    x, me_ex = get_me_exchange_array(sim_loc, ale, unmodeled)
                    me_ex = filter_exchange_fluxes(me_ex, threshold=.05)
                    make_exchange_met_plots(me_ex, x, unmodeled, ale,
                                            m_model=False, axis=exchange_ax)
                    exchange_ax.set_title(' & '.join([rxn_to_gene[i] for i in ale.split('-')]))

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
                         label='%s' % (suffix), **line_kwargs)

            if plot_kind == 'default':
                #m_gr = []
                #m_ex = defaultdict(list)
                #for frac in x:
                #    ko1s = ale.split('-')[0].split(':')
                #    ko2s = ale.split('-')[1].split(':')
                #    com_m = make_community_m_model(ijo, frac, ko1s, ko2s,
                #                                   scale_uptake=False)
                #    embed()
                #    m_sol = com_m.optimize()
                #    m_gr.append(m_sol.f)
                #    m_ex = get_metabolite_exchange(m_sol.x_dict, m_ex, frac)
#
                #if normalize_to_max:
                #    m_gr = np.array(m_gr) / np.array(m_gr).max()
#
                #m_ex = filter_exchange_fluxes(m_ex)
                #make_exchange_met_plots(m_ex, x, unmodeled, ale, m_model=True)
#
                #ax1.plot(x, m_gr, '--', label='M-model', marker='o',
                #         linewidth=.5)
#

                # grey region of max growth on plot
                me_max_loc = np.argmax(me_gr)
                ymin, ymax = ax1.get_ylim()
                ax1.fill_between([x[me_max_loc-1], x[me_max_loc], x[me_max_loc+1]],
                                 [me_gr[me_max_loc-1], me_gr[me_max_loc], me_gr[me_max_loc+1]],
                                 [ymin, ymin, ymin], facecolor='#B6C3C5', alpha=.3)
                ax1.plot([x[me_max_loc], x[me_max_loc]],
                         [ymin, me_gr[me_max_loc]], 'k--')
                # plotting changes ymin, set again
                ax1.set_ylim([ymin, ymax])

            if plot_kind != 'metabolite_limitation':
                lgd = axes1[1].legend(ncol=2, frameon=True, title=r'$\mathbf{\frac{K_{eff} \ Export \ \Delta hisD}{K_{eff} \ Export \ Partner \ Strain}}$')
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
                lgd = ax1.legend(sorted_handles, sorted_labels, ncol=1,
                                 frameon=True,
                                 title=r'$\mathbf{Metabolite \ Cross\mathrm{-}fed}$')
                if lgd:
                    plt.setp(lgd.get_title(), fontsize=12)

        if plot and subplot_axes is None:
            axes1[0].set_ylabel('Community Growth Rate')
            ax1.set_xlim([0, 1.])
            ax1.set_title(' & '.join([rxn_to_gene[i] for i in ale.split('-')]))
            fig1.tight_layout()
            fig1.savefig('%s/community_plots/%s/%s_%s.png' % (
                here, plot_kind, save_prefix, ale))

    if single_exchange_plot:
        exchange_axes[0].set_ylabel(r'Cross-feeding Flux ($\frac{mmol}{gDW_{shared}^{-1} \cdot hr^{-1}}$)')
        exchange_axes[1].set_xlabel('Fraction Strain 1')
        exchange_fig.tight_layout()
        exchange_fig.savefig('/home/sbrg-cjlloyd/Dropbox/%s_exchange.png' % plot_kind)


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
                      size=5, color="0", linewidth=0, ax=ax, orient='horizontal')
        ax.set_yticks([0, 1, 2, 3, 4])

        ax.set_yticklabels(["Metabolite restriction \n" + r"($\mathbf{Panel \ A}$)",
                            r"$\Delta$hisD" + " export less efficient \n " + r"($\mathbf{Panel \ B, \ solid}$)",
                            r"$\Delta$hisD" + " export more efficient \n " + r"($\mathbf{Panel \ B, \ dashed}$)",
                            "Experimentally Inferred \n " + r"($\mathbf{By \ mutation}$)",
                            "Experimentally Inferred \n " + r"($\mathbf{By \ coverage}$)"],
                           va='center', ha='center')
        ax.tick_params(axis='y', which='major', pad=90)
        ax.set_xlim([0, 1.])
        # ax.set_title(' & '.join([rxn_to_gene[i] for i in ale.split('-')]))
    axes[1].set_xlabel(r"$\Delta$hisD Fraction", fontsize=23,
                       color='black')

    if subplot_axes is None:
        fig.tight_layout()
        fig.savefig('%s/community_plots/%s/sweeps_box_plot.png' % (
            here, plot_kind))


if __name__ == '__main__':
    abundance_df = \
        pandas.read_excel('%s/relative_abundance/tables/characteristic_abundances_strain1.xlsx' % here,
                          index_col=0, header=None)

    plot_kinds = ['default', 'secretion_keff_sweep',
                  'metabolite_limitation', 'unmodeled_sweep']

    plot_kinds = ['default', 'metabolite_limitation', 'secretion_keff_sweep']
    normalize_to_max = False
    main_text_figure = True

    if main_text_figure:
        max_gr_dict = {}
        fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharey='row', sharex='col')
        for i, plot_kind in enumerate(plot_kinds):
            save_loc = '%s/community_plots/%s' % (here, plot_kind)
            if not os.path.exists(save_loc):
                os.mkdir(save_loc)

            sim_loc = expanduser('~/final_community_sims/%s' % plot_kind)
            make_computational_abundance_plots(abundance_df, sim_loc, plot_kind,
                                               normalize_to_max, max_gr_dict,
                                               subplot_axes=axes[i])
            print(max_gr_dict)

        make_sweeps_box_plots(max_gr_dict, abundance_df, plot_kind,
                              subplot_axes=axes[2])

        for i, ale in enumerate(['HISTD-CS', 'HISTD-DHORTS',
                                 'HISTD-GLUDy:GLUSy']):
            axes[0][i].set_title(' & '.join([rxn_to_gene[i] for
                                             i in ale.split('-')]),
                                 fontsize=25)

        # Add labels to each panel
        for i, t in enumerate(['A', 'B', 'C']):
            fig.text(0.085, 1-i/3.2, t, fontsize=35, fontweight='bold', va='top',
                     color='black')

        fig.tight_layout()
        fig.savefig('/home/sbrg-cjlloyd/Dropbox/Figure_8.png', dpi=250)

    else:
        plt.rcParams['axes.facecolor'] = 'w'
        max_gr_dict = {}
        plot_kind = 'default'
        save_loc = '%s/community_plots/%s' % (here, plot_kind)
        if not os.path.exists(save_loc):
            os.mkdir(save_loc)

        sim_loc = expanduser('~/final_community_sims/%s' % plot_kind)
        make_computational_abundance_plots(abundance_df, sim_loc, plot_kind,
                                           normalize_to_max, max_gr_dict,
                                           subplot_axes=None)
        print(max_gr_dict)
