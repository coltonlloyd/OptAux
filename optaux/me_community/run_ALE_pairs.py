from __future__ import print_function, division

from qminospy.me1 import ME_NLP1
import os
import numpy as np
import pandas as pd
import pickle
import cobra
import json

from optaux.resources.possible_uptake import get_possible_uptake
from optaux.me_community.me_model_community import (
    knock_out_reactions,
    change_unmodeled_protein_fraction, change_secretion_keff,
    scale_exchange_fluxes)

here = os.path.dirname(os.path.abspath(__file__))
resource_dir = here.replace('ME_community', 'resources')
ijo = cobra.io.load_json_model(os.path.join(resource_dir, 'iJO1366.json'))

# make IDs compatible with bigg ids
df = pd.read_csv(os.path.join(resource_dir, 'bigg_model_changes.csv'))
df = df[['old_reaction', 'new_reaction']].set_index('old_reaction')
bigg_id_change = {}
for key, value in df.itertuples():
    if type(key) is float:
        continue
    bigg_id_change[key] = value

# some exchange reactions are unique to ME-models, remove
skip_list = ['LI_c', 'pqq_e', 'cs_e', 'tl_c', 'cu_c', 'C10H8O5_c', 'C9H9O4_c',
             'NiFeCoCN2_c', 'RNase_m5', 'RNase_m16', 'RNase_m23']

with open(os.path.join(resource_dir, 'need_exporters.json'), 'r') as f:
    need_exporters = json.load(f)


def restrict_uptake(model, ko1, ko2, is_me_model=True):
    """If metabolite cannot restore growth to auxotroph,
    do not allow cross-feeding of it. Works for both M and ME-model"""

    for rxn in model.reactions.query('EX_'):
        if not rxn.id.startswith('EX_'):
            continue
        elif '_Shared' in rxn.id:
            continue
        if '_S1' in rxn.id:
            me_rxn = rxn.id.replace('_S1', '')
        elif '_S2' in rxn.id:
            me_rxn = rxn.id.replace('_S2', '')
        else:
            raise UserWarning('%s does not have _S1 or _S2 in id' % rxn.id)
        # make IDs compatible with bigg ids
        me_rxn = bigg_id_change.get(me_rxn, me_rxn)

        # some exchange reactions are unique to ME-models, remove
        if me_rxn.replace('EX_', '') in skip_list:
            continue
        ijo_rxn = ijo.reactions.get_by_id(me_rxn)

        aux_met = get_possible_uptake(ijo, ko1) if '_S1' in rxn.id \
            else get_possible_uptake(ijo, ko2)

        if me_rxn not in aux_met and ijo_rxn.lower_bound == 0.:
            # This function is ran before scaling exchange fluxes
            # (i.e. reactions are not split into forward and reverse yet)
            rxn.lower_bound = 0


def setup_simulation(me_model, ko1_list, ko2_list, fraction_strain_1,
                     unmodeled_protein_fractions=list(),
                     restrict_uptake_flag=False,
                     secretion_reaction_keffs=None, secretion_reactions=None,
                     scale_secretion=False, scale_uptake=False):

    print('With exporters for all mets without _c to _p export, but '
          'have _p to _c import')
    for m in need_exporters:
        for strain in ['S1', 'S2']:
            r = cobra.Reaction('%s_export_%s' % (m, strain))
            me_model.add_reaction(r)
            r.add_metabolites({'%s_c_%s' % (m, strain): -1,
                               '%s_p_%s' % (m, strain): 1})

    if restrict_uptake_flag:
        restrict_uptake(me_model, ko1_list, ko2_list)

    knock_out_reactions(me_model, ko1_list, ko2_list)
    # apply_fractional_abundance_growth_rate(me_model, fraction_strain_1)
    scale_exchange_fluxes(me_model, fraction_strain_1,
                          secretion=scale_secretion, uptake=scale_uptake)
    change_unmodeled_protein_fraction(
        me_model, pair_unmodeled_protein_fraction=unmodeled_protein_fractions)

    if secretion_reaction_keffs:
        print(secretion_reactions[1], secretion_reactions[0])
        print(secretion_reaction_keffs[0], secretion_reaction_keffs[1])
        # Change order of secretion reactions. Secretes the metabolite the other
        # strain requires
        change_secretion_keff(me_model, secretion_reactions[1],
                              secretion_reactions[0], secretion_reaction_keffs[0],
                              secretion_reaction_keffs[1])


def solve_community_model(me_model, hs, precision=1e-2, mumax=2.5):
    me_nlp = ME_NLP1(me_model, growth_key='mu')

    muopt, hs, xopt, cache = me_nlp.bisectmu(precision=precision, mumax=mumax,
                                             basis=hs)
    # Access the solution that is saved in the original minime object
    sol = me_model.solution
    sol.f = sol.x_dict['biomass_dilution_S1'] + sol.x_dict[
            'biomass_dilution_S2']
    return hs


def get_metabolic_flux(me_model):
    metabolic_flux = {}
    for r in ijo.reactions:
        for me_r in me_model.reactions.query(r.id):
            if '_Shared' in me_r.id:
                continue
            suffix = '_S1' if '_S1' in me_r.id else '_S2'
            sign = -1 if '_REV_' in me_r.id else 1
            if abs(me_r.x) > 1e-10:
                metabolic_flux[r.id + suffix] = sign * me_r.x
    return metabolic_flux


if __name__ == '__main__':
    ko1s = [['CS'], ['CS'], ['CS'], ['HISTD'], ['HISTD'], ['DHORTS']]
    ko2s = [['GLUDy', 'GLUSy'], ['HISTD'], ['DHORTS'],
            ['GLUDy', 'GLUSy'], ['DHORTS'], ['GLUDy', 'GLUSy']]

    unmodeled_fraction_list = [.4]
    fraction_strain1_list = np.linspace(.1, .9, 10)

    ko1s = [['HISTD'], ['HISTD'], ['HISTD']]
    ko2s = [['DHORTS'], ['GLUDy', 'GLUSy'], ['CS']]

   # ko1s = [['HISTD']]
   # ko2s = [['CS']]

    scale_secretion = True
    scale_uptake = True
    restrict_uptake_flag = False

    base_loc = os.path.join(here, 'no_abundance_based_unmodeled_scale_uptake_redo')
    for ko1_list, ko2_list in zip(ko1s, ko2s):
        if 'HISTD' not in ko1_list and 'HISTD' not in ko2_list:
            continue
        # initialize basis for solving
        hs = None

        pair_name = '-'.join(ko1_list) + ':' + '-'.join(ko2_list)
        print('Running simulations for %s pair' % pair_name)

        dir_name = '/'.join([base_loc, pair_name])
        if not os.path.isdir(dir_name):
            os.mkdir(dir_name)

        for q in unmodeled_fraction_list:
            fraction_dir = '%i_unmodeled_protein' % (q * 100.)
            output_dir = '/'.join([base_loc, pair_name, fraction_dir])
            print('Saving simulations in %s' % output_dir)

            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)

            sol_list, gr_list = [], []

            for fraction in fraction_strain1_list:
                with open(os.path.join(resource_dir,
                                       "iJL1678b_community.pickle"), 'rb') as f:
                    me = pickle.load(f)

                with open('./basis.pickle', 'rb') as f:
                    hs = pickle.load(f)
                setup_simulation(me, ko1_list, ko2_list, fraction,
                                 unmodeled_protein_fractions=[q, q],
                                 scale_secretion=scale_secretion,
                                 scale_uptake=scale_uptake,
                                 restrict_uptake_flag=restrict_uptake_flag)
                hs = solve_community_model(me, hs, precision=1e-4, mumax=2.5)
                with open('./basis.pickle', 'wb') as f:
                    pickle.dump(hs, f)
                with open(output_dir + '/%.2f_frac_strain1.pickle' %
                          fraction, 'wb') as f:
                    pickle.dump(me.solution, f)
                with open(output_dir + '/%.2f_frac_strain1_flux.json' %
                          fraction, 'w') as f:
                    json.dump(get_metabolic_flux(me), f)
