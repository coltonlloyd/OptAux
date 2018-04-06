from __future__ import print_function, division

import cobra
import cobra.test

all_kos = [['FUM', 'PPC', 'MALS'], ['FUM', 'PPC'], ['FUM', 'PPC', 'PSERT'], ['PSERT', 'PPC', 'SUCDi'],
           ['MALS', 'PPC'], ['ACONTa'], ['ICDHyr'], ['CS']]


def print_mets_for_auxs_fig_3():
    for kos in all_kos:
        ijo = cobra.test.create_test_model('ecoli')
        for ko in kos:
            ijo.reactions.get_by_id(ko).knock_out()
        for r in ['EX_fum_e', 'EX_succ_e', 'EX_glyclt_e',
                  'EX_mal__L_e', 'EX_cit_e', 'EX_akg_e', 'EX_asn__L_e']:
            substrate_rxn = ijo.reactions.get_by_id(r)
            substrate_rxn.lower_bound = - 10
            ijo.optimize()

            if ijo.solution.f > .01:
                print(kos, substrate_rxn.id, substrate_rxn.name)
            substrate_rxn.lower_bound = 0

if __name__ == '__main__':
    print_mets_for_auxs_fig_3()
