from __future__ import print_function, division

import json
import os

import pandas as pd
import cobra

from optaux.resources.possible_uptake import get_possible_uptake
from optaux.me_community.me_model_community import (
    knock_out_reactions,
    change_unmodeled_protein_fraction, change_secretion_keff,
    scale_exchange_fluxes)
from optaux import resources

here = os.path.dirname(os.path.abspath(__file__))
resource_dir = resources.__path__[0]
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
