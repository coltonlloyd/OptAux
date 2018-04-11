from __future__ import print_function, absolute_import, division

import pandas as pd

from optaux.helper_functions.characterize_auxotrophs import \
    get_auxotrophic_mets_per_ko
from optaux import resources

resource_dir = resources.__path__[0]
# Metabolites that can restore growth for each KO
ko_uptakes = {"['CS']": [u'EX_LalaDglu_e',
                         u'EX_LalaDgluMdap_e',
                         u'EX_LalaDgluMdapDala_e',
                         u'EX_LalaLglu_e',
                         u'EX_akg_e',
                         u'EX_arg__L_e',
                         u'EX_cit_e',
                         u'EX_fe3dcit_e',
                         u'EX_gln__L_e',
                         u'EX_glu__L_e',
                         u'EX_gthrd_e',
                         u'EX_orn_e',
                         u'EX_pro__L_e',
                         u'EX_progly_e'],
              "['DHORTS']": [u'EX_23ccmp_e',
                             u'EX_23cump_e',
                             u'EX_3cmp_e',
                             u'EX_3ump_e',
                             u'EX_cmp_e',
                             u'EX_csn_e',
                             u'EX_cytd_e',
                             u'EX_dcmp_e',
                             u'EX_dcyt_e',
                             u'EX_dump_e',
                             u'EX_duri_e',
                             u'EX_orot_e',
                             u'EX_uacgam_e',
                             u'EX_udpacgal_e',
                             u'EX_udpg_e',
                             u'EX_udpgal_e',
                             u'EX_udpglcur_e',
                             u'EX_ump_e',
                             u'EX_ura_e',
                             u'EX_uri_e'],
              "['GLUDy', 'GLUSy']": [u'EX_4abut_e',
                                     u'EX_LalaDglu_e',
                                     u'EX_LalaDgluMdap_e',
                                     u'EX_LalaDgluMdapDala_e',
                                     u'EX_LalaLglu_e',
                                     u'EX_agm_e',
                                     u'EX_ala__D_e',
                                     u'EX_ala__L_e',
                                     u'EX_alaala_e',
                                     u'EX_arg__L_e',
                                     u'EX_asn__L_e',
                                     u'EX_asp__L_e',
                                     u'EX_gln__L_e',
                                     u'EX_glu__L_e',
                                     u'EX_gthrd_e',
                                     u'EX_orn_e',
                                     u'EX_pro__L_e',
                                     u'EX_progly_e',
                                     u'EX_ptrc_e'],
              "['HISTD']": [u'EX_his__L_e'],
              "['ACGS']": ['EX_arg__L_e', 'EX_orn_e'],
              "['PPND']": ['EX_tyr__L_e', 'EX_tyrp_e'],
              "['PPNDH']": ['EX_phe__L_e'],
              "['DAPDC']": ['EX_frulys_e', 'EX_lys__L_e', 'EX_psclys_e'],
              "['HSST']": ['EX_met__L_e', 'EX_metsox__R__L_e',
                           'EX_metsox__S__L_e'],
              "['IGPDH', 'HISTP']": ['EX_his__L_e'],
              "['IPMD']": ['EX_leu__L_e'],
              "['PRAIi', 'IGPS']": ['EX_indole_e', 'EX_trp__L_e'],
              "['SERAT']": ['EX_cgly_e', 'EX_cys__L_e', 'EX_gthrd_e'],
              "['THRS', '4HTHRS']": ['EX_thr__L_e', 'EX_thrp_e']
              }


def get_possible_uptake(model, kos):
    # dict to correct difference in iJO/iJL1678-ME IDs
    full_map_df = pd.read_csv('%s/bigg_model_changes.csv' % resource_dir)
    map_df = full_map_df[['old_reaction', 'new_reaction']]
    map_dict = \
        map_df.set_index('new_reaction').dropna().to_dict()['old_reaction']

    if str(kos) in ko_uptakes:
        return ko_uptakes[str(kos)]
    elif str(list(reversed(kos))) in ko_uptakes:
        return ko_uptakes[str(list(reversed(kos)))]
    else:
        uptakes = get_auxotrophic_mets_per_ko(model, kos)
        return [map_dict.get(i, i) for i in uptakes]
