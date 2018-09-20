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
