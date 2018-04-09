import cobra
import pandas as pd
from optaux import resources
from optaux.core.characterize_auxotrophs import (
    get_auxotrophic_mets_per_ko, get_blocked_biomass, and_join_strings,
    to_string, gene_names_per_kos, genes_per_kos, id_to_name)
resource_dir = resources.__path__[0]

iJO = cobra.io.load_json_model('%s/iJO1366.json' % resource_dir)


if __name__ == '__main__':
    df = pd.DataFrame(columns=['Reaction Knockouts', 'Gene Knockouts',
                               'Gene Knockout Names',
                               'Blocked Biomass (BIGG)',
                               'Blocked Biomass (Names)',
                               'Auxotrophic Metabolites (BIGG)',
                               'Auxotrophic Metabolites (Names)'])
    for KOs in [['GLUDy', 'GLUSy'], ['CS'], ['HISTD'], ['DHORTS']]:
        auxotrophic_mets = get_auxotrophic_mets_per_ko(iJO, KOs)
        blocked_biomass = get_blocked_biomass(iJO, KOs)
        gene_names = and_join_strings(gene_names_per_kos(iJO, KOs))
        genes = and_join_strings(genes_per_kos(iJO, KOs))
        df_dict = {'Reaction Knockouts': and_join_strings(KOs),
                   'Gene Knockouts': genes,
                   'Gene Knockout Names': gene_names,
                   'Blocked Biomass (BIGG)': to_string(blocked_biomass),
                   'Blocked Biomass (Names)':
                       to_string(id_to_name(iJO, blocked_biomass)),
                   'Auxotrophic Metabolites (BIGG)':
                       to_string(auxotrophic_mets),
                   'Auxotrophic Metabolites (Names)':
                       to_string(id_to_name(iJO, auxotrophic_mets))
                   }
        df = df.append(pd.Series(df_dict), ignore_index=True)
        print(df)
    df.to_excel('experimental_auxotroph_characteristics.xls')
