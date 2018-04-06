import pandas as pd
from os.path import expanduser, dirname

from OptAux.core.characterize_auxotrophs import get_avg_flux_required
from OptAux import resources
import cobra
import numpy as np

xls = pd.ExcelFile(expanduser('~/Dropbox (UCSD SBRG)/'
                                '_OptAux_manuscript_and_scripts/Manuscript/'
                                'Supplements/'
                                'Supplementary Data File 1 - OptAux Solutions.xlsx'))

ijo = cobra.io.load_json_model(dirname(resources.__file__) + '/iJO1366.json')

ebc_kos = set()
mse_kos = set()
a_dict = {}
for sheet in xls.sheet_names:
    if '0_' not in sheet and '2_' not in sheet:
        continue

    df = xls.parse(sheet)

    df = df[df['Minimum Uptake at set_biomass'] < 0]
    for entry in df.index:
        if int(df.loc[entry, 'Number Auxotrophic Metabolites']) == 1:
            ebc_kos.add(df.loc[entry, 'Reaction Knockouts'])
        elif int(df.loc[entry, 'Number Auxotrophic Metabolites']) >= 5:
            mse_kos.add(df.loc[entry, 'Reaction Knockouts'])

        a_dict[df.loc[entry, 'Reaction Knockouts']] = df.loc[entry, 'Auxotrophic Metabolites (BIGG)']

ebc_avgs = set()
mse_avgs = set()

for kos in ebc_kos:
    aux_mets = a_dict[kos]

    ko_list = kos.split(' & ')
    avg = get_avg_flux_required(ijo, ko_list, [aux_mets])
    if avg:
        ebc_avgs.add(avg)

for kos in mse_kos:
    aux_mets = a_dict[kos]

    ko_list = kos.split(' & ')
    avg = get_avg_flux_required(ijo, ko_list, aux_mets.split(', '))
    if avg:
        mse_avgs.add(avg)

print(np.array(list(ebc_avgs)).mean())
print(np.array(list(mse_avgs)).mean())
