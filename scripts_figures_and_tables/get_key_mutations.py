import pandas as pd

from glob import glob

from optaux import resources


resource_dir = resources.__path__[0]
writer = pd.ExcelWriter('key_mutations.xls', engine='openpyxl')
for resequencing in glob('%s/resequencing_data/*Table.csv' % resource_dir):
    pair = '_'.join(resequencing.split('/')[-1].split(' ')[1:3])

    # load df
    df = pd.read_csv(resequencing)

    #df = df[df.columns.drop(list(df.filter(regex='I30 R1')))]
    #df = df[df.columns.drop(list(df.filter(regex='I40 R1')))]

    # make df of only mutation frequency values
    df_w_out_details = df[[i for i in df.columns if i.startswith('AUX')]]

    # only keep mutations that are found in at least one starting/ending clone
    df_key_mutations = \
        df_w_out_details[df_w_out_details > .9].dropna(how='all')

    # If mutation is called in one or more sample, it is deemed significant
    df_out = pd.DataFrame(columns=['Position', 'Mutation Type',
                                   'Sequence Change', 'Gene', 'Details'])

    for i in df_key_mutations.index:
        if df_w_out_details.loc[i].sum() >= 1:
            df_out = df_out.append(df.loc[i])

    df_out.to_excel(writer, sheet_name=pair)
writer.save()
