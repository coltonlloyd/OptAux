from __future__ import print_function, absolute_import, division

import re
import cobrame
from cobrame import mu
from sympy import Basic
import pickle
from os.path import abspath, dirname

from optaux import resources

here = dirname(abspath(__file__))
resource_dir = dirname(abspath(resources.__file__))


def make_binary_community_me(model, model_cons, savefile_name, ME_model=True):
    # Add 1 metabolite to model for each strain
    met_list = []
    for met in model.metabolites:
        met2 = met.copy()
        met.id += '_S1'
        met2.id += '_S2'
        met_list.append(met2)
    model.add_metabolites(met_list)
    model.repair()

    # Add 1 reaction to model for each strain
    new_rxn_list = []
    for rxn in model.reactions:
        rxn_type = type(rxn)
        new_rxn = rxn_type(rxn.id + '_S2')
        new_rxn.lower_bound = rxn.lower_bound
        new_rxn.upper_bound = rxn.upper_bound

        new_met_dict = {}
        for met, value in rxn.metabolites.items():
            new_met = model.metabolites.get_by_id(met.id.replace('_S1', '_S2'))
            new_met_dict[new_met] = value
        new_rxn.add_metabolites(new_met_dict)
        new_rxn_list.append(new_rxn)
        rxn.id += '_S1'
    model.add_reactions(new_rxn_list)
    model.repair()

    # Add exchanges into shared compartment for all exchanged metabolites
    media_list = []
    for rxn in model.reactions.query(re.compile('^EX_')):
        for met in rxn.metabolites:
            new_met = cobrame.Metabolite(
                met.id.replace('_S1', '_Shared').replace('_S2', '_Shared'))
            rxn.add_metabolites({new_met: 1})
        old_lower_bound = rxn.lower_bound
        rxn.lower_bound = -1000

        if '_S1' in rxn.id:
            new_rxn = cobrame.MEReaction(rxn.id.replace('_S1', '_Shared'))
            met = model.metabolites.get_by_id(new_rxn.id.replace('EX_', ''))
            new_rxn.add_metabolites({met: -1})
            new_rxn.upper_bound = 1000.
            new_rxn.lower_bound = 0.

            if rxn.lower_bound < 0:
                new_rxn.lower_bound = old_lower_bound * 2.
            model.add_reaction(new_rxn)

    for rxn in model_cons.reactions:
        new_rxn = model.reactions.get_by_id(rxn.id + '_S1')
        new_rxn_2 = model.reactions.get_by_id(rxn.id + '_S2')
        for met, value in rxn.metabolites.items():
            new_met = model.metabolites.get_by_id(met.id + '_S1')
            new_value = new_rxn.metabolites[new_met]
            if new_value != value:
                print('coeff different', rxn.id, met.id, new_value, value)

            new_met = model.metabolites.get_by_id(met.id + '_S2')
            new_value_2 = new_rxn_2.metabolites[new_met]
            if new_value_2 != value:
                print('coeff different', rxn.id, met.id, new_value_2, value)

            if not rxn.id.startswith('EX_') and (
                    new_rxn.lower_bound != rxn.lower_bound or
                    new_rxn_2.lower_bound != rxn.lower_bound):
                print('bounds different', rxn.id)
            if new_rxn.upper_bound != rxn.upper_bound or \
                    new_rxn_2.upper_bound != rxn.upper_bound:
                print('bounds different', rxn.id)

    if ME_model:
        # Add dummy objectives
        model.reactions.dummy_reaction_FWD_SPONT_S1.objective_coefficient = 1.
        model.reactions.dummy_reaction_FWD_SPONT_S2.objective_coefficient = 1.

        with open('%s/%s' % (resource_dir, savefile_name), 'wb') as f:
            pickle.dump(model, f)
    else:
        return model


def knock_out_reactions(model, ko1, ko2):
    for ko in ko1:
        for rxn in model.process_data.get_by_id(ko)._parent_reactions:
            model.reactions.get_by_id(rxn + '_S1').knock_out()
    for ko in ko2:
        for rxn in model.process_data.get_by_id(ko)._parent_reactions:
            model.reactions.get_by_id(rxn + '_S2').knock_out()


def scale_exchange_fluxes(model, fraction_S1, secretion=True, uptake=True):
    if secretion and uptake:
        print('Scaling Secretion and Uptake Flux')
    elif secretion:
        print('Scaling Secretion')
    elif uptake:
        print('Scaling Uptake')
    else:
        print('Not scaling exchange fluxes')

    for rxn in model.reactions.query('EX_'):
        # filter out a few complexes with 'COMPLEX_mod' in ID
        if not rxn.id.startswith('EX_') or 'Shared' in rxn.id:
            continue

        check_rxn_id = rxn.id.replace('_S1', '_Shared').replace('_S2',
                                                                '_Shared')
        check_rxn = model.reactions.get_by_id(check_rxn_id)
        if check_rxn.lower_bound < 0 and check_rxn.id != 'EX_glc__D_e_Shared':
            continue

        # split exchanges into forward and reverse reactions
        if rxn.lower_bound < 0:
            rev_rxn = cobrame.MEReaction(rxn.id + '_reverse')
            model.add_reaction(rev_rxn)
            rev_rxn.add_metabolites(rxn._metabolites)
            rev_rxn.lower_bound = -1000
            rev_rxn.upper_bound = 0
            rxn.lower_bound = 0

        if '_S1' not in rxn.id and '_S2' not in rxn.id:
            continue

        met_id = rxn.id.replace('EX_', '')
        strain = '_S1' if '_S1' in rxn.id else '_S2'
        met_shared_id = met_id.replace(strain, '_Shared')

        # Determine how to handle reaction stoichiometry change
        frac_term = fraction_S1 if strain is '_S1' else (1 - fraction_S1)

        if secretion:
            # met -> x * met_shared
            rxn.add_metabolites(
                {model.metabolites.get_by_id(met_shared_id): frac_term},
                combine=False)
            if rxn.id in ['EX_glc__D_e_S1', 'EX_glc__D_e_S2']:
                rxn.upper_bound = 0
                print(rxn.reaction, rxn.upper_bound)

            # met <- x * met_shared
            rxn_rev_id = rxn.id + '_reverse'
            if rxn_rev_id in model.reactions:
                rxn_rev = model.reactions.get_by_id(rxn_rev_id)
                rxn_rev.add_metabolites(
                    {model.metabolites.get_by_id(met_shared_id): (frac_term)},
                    combine=False)


def change_unmodeled_protein_fraction(model, fraction_S1=0,
                                      pair_unmodeled_protein_fraction=[],
                                      max_unmodeled=.9, min_unmodeled=None):
    """
    max_unmodeled = highest unmodeled protein fraction that ME-model can grow
    min_unmodeled = unmodeled protein fraction value when fraction of a strain
    """

    if pair_unmodeled_protein_fraction:
        print('User defined unmodeled protein fraction')
        value_1, value_2 = pair_unmodeled_protein_fraction
    else:
        print('Scaling unmodled protein fraction based on strain abundance and'
              'min-max defined values')
        # Low strain abundance -> low growth rate and high unmodeled
        # protein fraction
        fraction_1 = 1 - fraction_S1
        fraction_2 = fraction_S1
        value_1 = fraction_1 * (max_unmodeled - min_unmodeled) + min_unmodeled
        value_2 = fraction_2 * (max_unmodeled - min_unmodeled) + min_unmodeled

    amount_1 = value_1 / (1 - value_1)
    amount_2 = value_2 / (1 - value_2)

    for strain in ['S1', 'S2']:
        rxn = model.reactions.get_by_id('protein_biomass_to_biomass_%s' % strain)

        dummy = model.metabolites.get_by_id('unmodeled_protein_biomass_%s' % strain)
        biomass = model.metabolites.get_by_id('biomass_%s' % strain)

        if strain == 'S1':
            rxn.add_metabolites({dummy: -amount_1, biomass: amount_1 + 1},
                                combine=False)

        if strain == 'S2':
            rxn.add_metabolites({dummy: -amount_2, biomass: amount_2 + 1},
                                combine=False)


def change_secretion_keff(model, reaction1, reaction2, keff1, keff2):
    for r in model.reactions.query(reaction1 + '_'):
        if '_S2' in r.id:
            continue

        # Assumes form "[reaction_id]_[FWDorREV]_[enzyme_id]_[S1orS2]
        if '_FWD_' in r.id:
            enzyme = r.id.split('_FWD_')[1]
        elif '_REV_' in r.id:
            enzyme = r.id.split('_REV_')[1]
        else:
            raise UserWarning('Bad')

        r.add_metabolites({enzyme: -cobrame.mu / keff1 / 3600.},
                          combine=False)

    for r in model.reactions.query(reaction2 + '_'):
        if '_S1' in r.id:
            continue

        # Assumes form "[reaction_id]_[FWDorREV]_[enzyme_id]_[S1orS2]
        if '_FWD_' in r.id:
            enzyme = r.id.split('_FWD_')[1]
        elif '_REV_' in r.id:
            enzyme = r.id.split('_REV_')[1]
        else:
            raise UserWarning('Bad')

        r.add_metabolites({enzyme: -cobrame.mu / keff2 / 3600.},
                          combine=False)
