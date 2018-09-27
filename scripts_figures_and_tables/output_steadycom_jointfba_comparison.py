from optaux.resources.update_mpl_rcparams import update_rcparams
import numpy as np
from optaux.me_community import steadycom
import cobra
from matplotlib import pyplot as plt
import os
from optaux import resources
from optaux.resources.possible_uptake import ko_uptakes, get_possible_uptake
from optaux.me_community.me_model_community import make_binary_community_me, scale_exchange_fluxes


rxn_to_gene = {'CS': r'$\Delta \mathit{gltA} \Delta \mathit{prpC}$',
               'HISTD': r'$\Delta \mathit{hisD}$',
               'GLUDy:GLUSy': r'$\Delta \mathit{gdhA} \Delta \mathit{gltB}$',
               'DHORTS': r'$\Delta \mathit{pyrC}$'}

resource_dir = resources.__path__[0]

update_rcparams()
plt.rcParams['axes.facecolor'] = 'w'


fig, axes = plt.subplots(1, 3, figsize=(15, 5))

ijo = cobra.io.load_json_model('%s/iJO1366.json' % resource_dir)
r = cobra.Reaction('tranport_orot')
ijo.add_reaction(r)
r.add_metabolites({'orot_c': -1, 'orot_e': 1})

ko1s = [['HISTD']] * 3
ko2s = [['CS'], ['DHORTS'], ['GLUDy', 'GLUSy']]
i = 0
for ko1, ko2 in zip(ko1s, ko2s):
    model = steadycom.run(ijo.copy(), ko1, ko2,
                          restrict_uptake_to_one_metabolite=False)
    model.objective = 'community_biomass_variable'
    steadycom.solve(model, solver='gurobi')
    fva = steadycom.run_fva_for_abundances(model, .9,
                                           # Round to prevent numerical errors
                                           round(model.solution.f, 4),
                                           solver='gurobi')
    y = []
    x = []
    for frac in reversed(sorted(fva)):
        abundance_dict = fva[frac]['abundance_strain_0']
        y = [frac] + y + [frac]
        x = [abundance_dict['minimum']] + x + [abundance_dict['maximum']]

    axes[i].plot(x, y, linewidth=5, label='SteadyCom')
    i += 1

# JointFBA community M-modeling
for i, pair in enumerate([['CS'], ['DHORTS'], ['GLUDy', 'GLUSy']]):
    s1s = [0]
    grs = [0]

    for s1 in np.linspace(.05, .95, 10):
        com = make_binary_community_me(ijo.copy(), ijo.copy(), '',
                                       ME_model=False)

        com.reactions.HISTD_S1.knock_out()
        for r in pair:
            com.reactions.get_by_id(r + '_S2').knock_out()
        biomass = com.reactions.BIOMASS_Ec_iJO1366_core_53p95M_S1 + com.reactions.BIOMASS_Ec_iJO1366_core_53p95M_S2
        biomass.id = 'biomass'
        com.add_reaction(biomass)
        com.objective = biomass
        s1s.append(s1)
        scale_exchange_fluxes(com, s1)
        for r in com.reactions:
            if '_reverse' in r.id:
                r.id = r.id.replace('_reverse', '_back')
        com.repair()
        com.reactions.EX_glc__D_e_Shared.lower_bound = -5
        com.reactions.EX_glc__D_e_S1.upper_bound = 0
        com.reactions.EX_glc__D_e_S2.upper_bound = 0
        for rxn in com.reactions.query('EX_'):
            if not rxn.id.startswith('EX_') or rxn.id.endswith(
                    '_S1') or rxn.id.endswith('_S2'):
                continue
            elif '_Shared' in rxn.id:
                continue
            if '_S1_back' in rxn.id:
                me_rxn = rxn.id.replace('_S1_back', '')
            elif '_S2_back' in rxn.id:
                me_rxn = rxn.id.replace('_S2_back', '')
            else:
                raise UserWarning('%s does not have _S1 or _S2 in id' % rxn.id)

            ijo_rxn = ijo.reactions.get_by_id(me_rxn)

            aux_met = 'EX_his__L_e' if '_S1' in rxn.id else ko_uptakes[
                str(pair)]

            if me_rxn not in aux_met and ijo_rxn.lower_bound == 0.:
                # This function is ran before scaling exchange fluxes
                # (i.e. reactions are not split into forward and reverse yet)
                rxn.lower_bound = 0

        sol = com.optimize()
        grs.append(sol.f)

    s1s.append(1)
    grs.append(0)
    axes[i].plot(s1s, np.array(grs) / np.array(grs).max(), '--',
                 label='M-model', linewidth=5)
    axes[i].set_title(r'%s & %s' % (rxn_to_gene['HISTD'],
                                    rxn_to_gene[':'.join(pair)]))
    # plt.show()
    axes[i].legend()
    axes[i].set_xlabel('Fraction Strain 1')
    axes[i].set_ylim((.9, 1.01))
axes[0].set_ylabel(r'Community Growth Rate (Normalized)')
fig.savefig('steadycom_m_model_comparison.png')