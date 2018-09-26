import pandas as pd
import cobra
from optaux.helper_functions import characterize_auxotrophs
from optaux import resources
import glob
import pandas as pd
from collections import Counter


resource_dir = resources.__path__[0]

ijo = cobra.io.load_json_model(resource_dir + '/iJO1366.json')

df = pd.DataFrame()
for fi in glob.glob('./optaux*.xls'):
    temp_df = pd.read_excel(fi)
    df = df.append(temp_df, ignore_index=True)


def output_mse_summary_df():
    i = 0
    df_out = pd.DataFrame()
    count = Counter()
    for a in df[df['Minimum Uptake at set_biomass'] < 0].groupby(
            by='Auxotrophic Metabolites (BIGG)'):
        kos = set([i for i in a[1]['Reaction Knockouts'].values])

        for ko_list in kos:
            i += 1
            subsystems = set()
            num_auxs = a[1]['Number Auxotrophic Metabolites'].values[0]
            if num_auxs < 5:
                continue
            for ko in ko_list.split(' & '):
                if ijo.reactions.get_by_id(ko).subsystem:
                    subsystems.add(ijo.reactions.get_by_id(ko).subsystem)
                else:
                    print(ijo.reactions.get_by_id(ko).subsystem)
            count.update(subsystems)
            sub_list = sorted(list(subsystems))
            print(ko_list, sub_list)
            df_out.loc[i, "Subsystem 1"] = sub_list[0]
            if len(sub_list) > 1:
                df_out.loc[i, "Subsystem 2"] = sub_list[1]
            if len(sub_list) > 2:
                df_out.loc[i, "Subsystem 3"] = sub_list[2]
            df_out.loc[i, "KOs"] = ko_list
            df_out.loc[i, 'Number Auxotrophic Metabolites'] = num_auxs
            df_out.loc[
                i, 'Avgerage Uptake Required Flux'] = characterize_auxotrophs.get_avg_flux_required(
                ijo, ko_list.split(' & '), a[0].split(', '))

    df_out.sort_values(['Subsystem 1', 'Subsystem 2',
                        'Subsystem 3']).to_excel('MSE_summary.xls')


def print_unique_and_total_designs():
    i = 0
    total = 0
    unique_total = 0
    print('---------------------All Designs-------------------')
    for a in df[df['Minimum Uptake at set_biomass'] < 0].groupby(
            by='Auxotrophic Metabolites (BIGG)'):
        kos = set()
        unique = set()
        i += 1
        kos.update([i for i in a[1]['Reaction Knockouts'].values])
        unique.update([str(i) for i in a[1]['Auxotrophic Metabolites (BIGG)']])
        total += len(kos)
        unique_total += len(unique)

    print('Number of unique auxotrophic designs:', i)
    print('Number of designs total:', total)

    print('------------------EBC Designs---------------------')
    i = 0
    total = 0

    for a in df[df['Minimum Uptake at set_biomass'] < 0].groupby(
            by='Auxotrophic Metabolites (BIGG)'):
        if a[1]['Number Auxotrophic Metabolites'].values[0] == 1:
            kos = set()
            unique = set()
            i += 1
            kos.update([i for i in a[1]['Reaction Knockouts'].values])
            total += len(kos)

    print('Number of unique auxotrophic designs:', i)
    print('Number of designs total:', total)

    print('------------------EBC Designs (semi-specific)---------------------')
    i = 0
    total = 0

    for a in df[df['Minimum Uptake at set_biomass'] < 0].groupby(
            by='Auxotrophic Metabolites (BIGG)'):
        if 1 < a[1]['Number Auxotrophic Metabolites'].values[0] < 5:
            kos = set()
            unique = set()
            i += 1
            kos.update([i for i in a[1]['Reaction Knockouts'].values])
            total += len(kos)

    print('Number of unique auxotrophic designs:', i)
    print('Number of designs total:', total)
    print('------------------MSE Designs ---------------------')
    i = 0
    total = 0

    for a in df[df['Minimum Uptake at set_biomass'] < 0].groupby(
            by='Auxotrophic Metabolites (BIGG)'):
        if a[1]['Number Auxotrophic Metabolites'].values[0] >= 5:
            kos = set()
            unique = set()
            i += 1
            kos.update([i for i in a[1]['Reaction Knockouts'].values])
            total += len(kos)

    print('Number of unique auxotrophic designs:', i)
    print('Number of designs total:', total)


def output_table_of_simplest_ebc_design():
    '''Prints desings with smallest number of KOs. Used to populate table
    in S1 appendix'''

    def keep_only_shortest_knockout(rxns):
        shortest_length = 5
        rxn_list = [r.split(' & ') for r in sorted(rxns, key=len)]
        for r in rxn_list:
            if len(r) < shortest_length:
                shortest_length = len(r)
        for r in list(rxn_list):
            if len(r) > shortest_length:
                rxn_list.remove(r)
        return set([' & '.join(rs) for rs in rxn_list])

    i = 0
    out_df = pd.DataFrame()
    for a in df[df['Minimum Uptake at set_biomass'] < 0].groupby(
            by='Auxotrophic Metabolites (BIGG)'):
        if a[1]['Number Auxotrophic Metabolites'].values[0] != 1:
            continue
        kos = set()
        rxns = set()
        genes = set()

        i += 1
        kos.update([i for i in a[1]['Gene Knockout Names'].values])
        genes.update([i for i in a[1]['Gene Knockouts'].values])
        rxns.update(
            [i for i in sorted(a[1]['Reaction Knockouts'].values, key=len)])
        rxns = keep_only_shortest_knockout(rxns.copy())
        met_name = ijo.reactions.get_by_id(a[0]).name.replace(' exchange', '')
        out_df.loc[met_name, 'Reactions'] = ''
        for rxn_set in rxns:
            for r in rxn_set.split(' & '):
                out_df.loc[met_name, 'Reactions'] += '%s (%s)' % (
                r, ijo.reactions.get_by_id(r).gene_name_reaction_rule)
                out_df.loc[met_name, 'Reactions'] += ' & '
            out_df.loc[met_name, 'Reactions'] = out_df.loc[
                met_name, 'Reactions'].rstrip(' & ')
            out_df.loc[met_name, 'Reactions'] += ' | '
        out_df.loc[met_name, 'Reactions'] = out_df.loc[
            met_name, 'Reactions'].rstrip(' | ')

    out_df.to_csv('ebc_simplest_designs.csv')


if __name__ == '__main__':

    output_table_of_simplest_ebc_design()
    print_unique_and_total_designs()
    output_mse_summary_df()