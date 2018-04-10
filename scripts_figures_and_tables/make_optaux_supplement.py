from __future__ import print_function, absolute_import, division

from datetime import date

import pandas as pd
import cobra

from optaux.core.run_optaux import run_optaux
from optaux import resources
from optaux.helper_functions.characterize_auxotrophs import (
    get_auxotrophic_mets_per_ko, get_blocked_biomass, and_join_strings,
    to_string, gene_names_per_kos, genes_per_kos, id_to_name)


if __name__ == '__main__':
    resource_dir = resources.__path__[0]
    model_loc = resource_dir + 'iJO1366.json'
    false_positive_loc = resource_dir + 'OrthKOEssential.csv'

    media_list = []

    lb_essential_list = \
        pd.read_csv(resource_dir + 'LB_essential_ecocyc_iJO1366.csv',
                    index_col=0).T.values[0]
    # Load model for OptAux
    cons_model = cobra.io.load_json_model(model_loc)
    writer = pd.ExcelWriter('supplement_1_optaux_solutions.xls')

    # Load known model false positives. Exclude these from possible knockouts

    false_positives = pd.read_csv(false_positive_loc,
                                  index_col=0).index.values

    for trace in [2]:
        df = pd.DataFrame(columns=['Target Metabolite',
                                   'Reaction Knockouts',
                                   'Gene Knockouts',
                                   'Gene Knockout Names',
                                   'Oxygen Uptake',
                                   'Number of Knockouts',
                                   'Trace Metabolite Threshold',
                                   'Growth Rate used for set_biomass',
                                   'Minimum Uptake at set_biomass',
                                   'Maximum Growth Rate w/ Knockouts',
                                   'Minimum Uptake at Max Growth Rate w/ Knockouts',
                                   'Time (s)', 'Date', 'Solver',
                                   'Solve Status',
                                   'Blocked Biomass (BIGG)',
                                   'Blocked Biomass (Names)',
                                   'Number Auxotrophic Metabolites',
                                   'Auxotrophic Metabolites (BIGG)',
                                   'Auxotrophic Metabolites (Names)'
                                   ])
        for num in [1, 2, 3]:
            for r in cons_model.reactions.query('EX_'):
                print(r.id)

                # Load model for OptAux
                model = cobra.io.load_json_model(model_loc)
                met = r.id.replace('EX_', '')

                # Fraction of maximum growth rate for optaux sim
                frac = .01
                if 'C' not in model.metabolites.get_by_id(met).elements:
                    continue
                try:
                    ans = run_optaux(model, met, frac, num,
                                     lb_essential_list=lb_essential_list,
                                     media_list=[], solver='gurobi',
                                     exclude_reactions=false_positives,
                                     trace_metabolite_threshold=trace)
                except:
                    print('Error in %i_%s' % (num, r.id))
                    continue
                if not ans:
                    continue

                # Add characterization of OptAux knockout solutions
                KOs = ans['Reaction Knockouts']
                auxotrophic_mets = get_auxotrophic_mets_per_ko(cons_model, KOs)
                blocked_biomass = get_blocked_biomass(cons_model, KOs)
                gene_names = and_join_strings(gene_names_per_kos(cons_model,
                                                                 KOs))
                genes = and_join_strings(genes_per_kos(cons_model, KOs))
                ans.update({'Reaction Knockouts': and_join_strings(KOs),
                            'Gene Knockouts': genes,
                            'Gene Knockout Names': gene_names,
                            'Blocked Biomass (BIGG)': to_string(
                                blocked_biomass),
                            'Blocked Biomass (Names)':
                                to_string(id_to_name(cons_model,
                                                     blocked_biomass)),
                            'Number Auxotrophic Metabolites':
                                len(auxotrophic_mets),
                            'Auxotrophic Metabolites (BIGG)':
                                to_string(auxotrophic_mets),
                            'Auxotrophic Metabolites (Names)':
                                to_string(id_to_name(cons_model,
                                                     auxotrophic_mets))
                            })

                ans['Target Metabolite'] = met
                ans['Date'] = str(date.today())

                df = df.append(ans, ignore_index=True)
                df.to_excel('optaux_intermediate_trace_%.2f.xls' % trace)
        df.to_excel(excel_writer=writer, sheet_name='%.2f Trace' % trace)
    writer.save()
