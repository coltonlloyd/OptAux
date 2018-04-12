from __future__ import print_function, absolute_import, division

from glob import glob
from collections import defaultdict
import json
import os
from cycler import cycler

import numpy as np
from matplotlib import pyplot as plt
from optaux.resources.update_mpl_rcparams import update_rcparams

update_rcparams()
colors = ['#0099E6', '#F23814', '#DFB78B', '#9E9E9E', '#12255A']
colors.extend(reversed(colors[:-1]))
plt.rcParams['axes.prop_cycle'] = cycler(color=colors)


skip_list = ['h2o', 'etoh', 'o2', 'mg2', 'nh4', 'so4', 'h', 'ac', 'for',
             'pi', 'lac__D', 'co2', 'glc__D', 'fe2', 'fe3', 'fum']

here = os.path.dirname(os.path.abspath(__file__))

rxn_to_gene = {'CS': r'$\Delta \mathit{gltA} \Delta \mathit{prpC}$',
               'HISTD': r'$\Delta \mathit{hisD}$',
               'GLUDy:GLUSy': r'$\Delta \mathit{gdhA} \Delta \mathit{gltB}$',
               'DHORTS': r'$\Delta \mathit{pyrC}$'}


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


def filter_exchange_fluxes(exchange_fluxes, sim_kind, threshold=.2):
    output_fluxes = {}
    for r, flux in exchange_fluxes.items():

        flux = np.array(flux)

        if len(flux) == 0 or abs(flux).sum() <= threshold:
            continue
        if sim_kind == 'default' and 'his__L' in r:
            continue
        output_fluxes[r] = flux

    return output_fluxes


def get_me_exchange_array(sim_loc, ale, unmodeled):
    me_ex = defaultdict(list)
    x = []

    simulations = '/'.join([sim_loc, ale, unmodeled, '*.json'])
    sims = glob(simulations)
    if len(sims) == 0:
        return None, None

    for sim in sorted(sims):
        # fraction of strain 1
        frac = (float(sim.split('/')[-1].split('_')[0]))
        if int(100. * frac) % 5 != 0 or int(100. * frac) % 10 == 0:
            continue
        x.append(frac)
        with open(sim, 'r') as f:
            x_dict = json.load(f)

        me_ex = get_metabolite_exchange(x_dict, me_ex, frac)

    return x, me_ex


def add_met_exchange_to_axis(ex, x, sim_kind, ale, ax):
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
        if sim_kind != 'default':
            color = colors[0] if met == 'his__L' else colors[1]
            ax.plot(x, flux, line_type, linewidth=5,
                    label=met, color=color)
        else:
            ax.plot(x, flux, line_type, linewidth=5, label=met)
        i += 1

    if i > 2:
        ncol = 2
    else:
        ncol = 1
    ax.legend(ncol=ncol, fontsize=15)
    ax.set_title(ale)
    ax.set_xlim([.05, .95])


def plot_metabolit_exchange(sim_loc, ales, unmodeled, sim_kind):
    fig, axes = plt.subplots(1, 3, sharex='row', sharey='col', figsize=(15, 4))
    for i, ale in enumerate(ales):
        ax = axes[i]
        x, me_ex = get_me_exchange_array(sim_loc, ale, unmodeled)
        me_ex = filter_exchange_fluxes(me_ex, sim_kind, threshold=.02)
        add_met_exchange_to_axis(me_ex, x, sim_kind, ale, ax)
        ax.set_title(
            ' & '.join([rxn_to_gene[i] for i in ale.split('-')]))

    axes[0].set_ylabel(
        r'Cross-feeding Flux ($\frac{mmol}{gDW_{shared}^{-1} \cdot hr^{-1}}$)')
    axes[1].set_xlabel(r'Fraction $\Delta \mathit{hisD}$')
    fig.tight_layout()
    fig.savefig('%s/community_plots/exchange_flux_plots/%s_exchange.png' %
                (here, sim_kind))


if __name__ == '__main__':
    plt.rcParams['axes.facecolor'] = 'w'
    ales = ['HISTD-CS', 'HISTD-DHORTS', 'HISTD-GLUDy:GLUSy']
    sim_location = '%s/community_me_sims/' % here
    save_location = '%s/community_plots/' % here

    sim_kind = 'secretion_keff_sweep'
    sim_loc = '/'.join([sim_location, sim_kind])
    # plot with default values
    unmodeled = '100.00_100.00_0_secretion_multiplier'
    plot_metabolit_exchange(sim_loc, ales, unmodeled, sim_kind)

    sim_kind = 'default'
    sim_loc = '/'.join([sim_location, sim_kind])
    unmodeled = '75.00_75.00_unmodeled_protein'
    plot_metabolit_exchange(sim_loc, ales, unmodeled, sim_kind)
