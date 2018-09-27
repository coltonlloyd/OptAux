import json
from matplotlib import pyplot as plt
import glob
from optaux.resources.update_mpl_rcparams import update_rcparams
import numpy as np
import cobra
from os.path import dirname, abspath
from matplotlib import pyplot as plt
from optaux import resources
from optaux.resources.possible_uptake import ko_uptakes
from optaux.me_community.me_model_community import make_binary_community_me, scale_exchange_fluxes

here = dirname(abspath(__file__))
resource_dir = dirname(abspath(resources.__file__))
sim_dir = '%s/community_sims_output_default_keffs/' % here
save_location = '%s/community_plots/' % here

update_rcparams()
plt.rcParams['axes.facecolor'] = 'w'


# ############## Plot glucose limited ME-model simulations ####################

pairs = ['HISTD-CS', 'HISTD-DHORTS', 'HISTD-GLUDy:GLUSy']
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for i, pair in enumerate(pairs):
    fracs = [0]
    grs = [0]
    for fi in sorted(glob.glob('%s/glucose_limited/%s/*/*_sol.json' %
                               (sim_dir, pair))):
        frac = float(fi.split('/')[-1].split('_')[0])
        with open(fi, 'r') as f:
            sol = json.load(f)
        gr = sol['biomass_dilution_S1']
        fracs.append(frac)
        grs.append(gr)
    grs.append(0)
    fracs.append(1)
    axes[i].plot(fracs, np.array(grs) / np.array(grs).max(),
                 label='ME-model (Substrate Limited)', linewidth=3,
                 marker='o')

    axes[i].set_title(pair)

# ############## Run and plot M-model simulations #############################

ijo = cobra.io.load_json_model('%s/iJO1366.json' % resource_dir)

# Add orotate transport for DHORTS knockout
r = cobra.Reaction('tranport_orot')
ijo.add_reaction(r)
r.add_metabolites({'orot_c': -1, 'orot_e': 1})

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
                 label='M-model', linewidth=3, marker='o')
    # axes[i].set_title(str(pair))
    # plt.show()
    axes[i].legend()
    axes[i].set_xlabel(r'Fraction Strain \Delta hisD}')
    axes[i].set_ylim((.9, 1.01))
    axes[i].set_xlim((0, 1))
axes[0].set_ylabel('Computed Community \n Growth Rate (Normalized)')
fig.savefig('%s/Figure_7_me_m_comparison.png' % save_location)
