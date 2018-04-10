from __future__ import print_function

from os.path import dirname, abspath
import pandas as pd
import numpy as np
import glob

import pysam
from Bio import SeqIO

from optaux import resources

resource_dir = resources.__path__[0]
here = dirname(abspath(__file__))

strain_to_characteristic = \
    {'hisD': ['2,400,257', '3,178,254'],  # lrhA/alaA, yqiC
     'gltB': ['3,178,277', '4,546,084'],  # yqiC, yjiC
     'gltA': ['2,719,426'],  # kgtP/rrfG
     'pyrC': ['444,779']  # cyoB
     }

pair_to_kos = {'hisD_pyrC': [['hisD'], ['pyrC']],
               'hisD_gltA': [['hisD'], ['gltA', 'prpC']],
               'hisD_gltB': [['hisD'], ['gltB', 'gdhA']]}

gene_to_reaction_ko_abbrev = {'hisD': 'HISTD',
                              'gltA': 'CS',
                              'pyrC': 'DHORTS',
                              'gltB': 'GLUDy'}


def _fill_missing_df_entries(df, pair):
    if pair != 'hisDpyrC':
        return df
    # lrhA, alaA (hisD characteristic mutation)
    df.loc['2,400,257', 'AUX hisD pyrC A3 F26 I0 R1'] = .144  # likely same mutation

    # cyoB (pyrC)
    df.loc['444,779', 'AUX hisD pyrC A3 F17 I0 R1'] = .905
    df.loc['444,779', 'AUX hisD pyrC A3 F26 I0 R2'] = .961
    df.loc['444,779', 'AUX hisD pyrC A4 F11 I0 R1'] = .836
    df.loc['444,779', 'AUX hisD pyrC A4 F11 I0 R2'] = .857
    df.loc['444,779', 'AUX hisD pyrC A4 F17 I0 R1'] = .892

    # yqiC (hisD characteristic mutation)
    df.loc['3,178,254', 'AUX hisD pyrC A2 F10 I0 R1'] = .292
    df.loc['3,178,254', 'AUX hisD pyrC A3 F26 I0 R1'] = .132
    df.loc['3,178,254', 'AUX hisD pyrC A3 F26 I0 R2'] = 0.  # not found
    df.loc['3,178,254', 'AUX hisD pyrC A4 F11 I0 R1'] = 0.  # not found
    df.loc['3,178,254', 'AUX hisD pyrC A4 F11 I0 R2'] = .167
    df.loc['3,178,254', 'AUX hisD pyrC A4 F17 I0 R1'] = .111
    df.loc['3,178,254', 'AUX hisD pyrC A4 F24 I0 R1'] = 0.  # not found
    df.loc['3,178,254', 'AUX hisD pyrC A4 F24 I0 R2'] = .088

    # envZ (shouldn't be 100%, breseq rounds values to 100 in some cases)
    df.loc['3,529,197', 'AUX hisD pyrC A4 F11 I0 R2'] = .719
    df.loc['3,529,197', 'AUX hisD pyrC A4 F17 I0 R1'] = .915
    df.loc['3,529,197', 'AUX hisD pyrC A4 F24 I0 R2'] = .913

    return df


def _add_ale_numbers_to_df(df):
    df = df.T
    for i, exp in enumerate(df.T):
        df.loc[exp, 'Ale'] = int(exp.split(' ')[3].replace('A', ''))
        df.loc[exp, 'Flask'] = int(exp.split(' ')[4].replace('F', ''))
        df.loc[exp, 'Isolate'] = int(exp.split(' ')[5].replace('I', ''))
    df.set_index(['Ale', 'Flask', 'Isolate'], inplace=True)
    df = df.T.sort_index(axis=1)
    return df


def _append_rows_to_abundance_df(out_df, strain1_name, strain2_name,
                                 strain1_abundance, strain2_abundance):

    # get reaction name from gene KO
    strain_1_ko = gene_to_reaction_ko_abbrev[strain1_name]
    strain_2_ko = gene_to_reaction_ko_abbrev[strain2_name]

    # normalize abundance to 100%
    strain_1 = strain1_abundance / (strain1_abundance + strain2_abundance)
    strain_2 = strain2_abundance / (strain1_abundance + strain2_abundance)

    out_df = \
        out_df.append(pd.Series(strain_1, name='__'.join([strain_1_ko,
                                                          strain_2_ko])))
    out_df = \
        out_df.append(pd.Series(strain_2, name='__'.join([strain_2_ko,
                                                          strain_1_ko])))
    return out_df


def _get_abundance_statistics(df, strain_2_name, stats_df):
    df_new = df.T.reset_index()
    for ale in set(df_new['Ale']):
        filtered_df = df_new[df_new['Ale'] == ale]
        analysis_series = filtered_df['hisD'] / (filtered_df['hisD'] +
                                                 filtered_df[strain_2_name])
        stats_df.loc[ale, 'mean'] = analysis_series.mean()
        stats_df.loc[ale, 'Stdev'] = analysis_series.std()

        print(ale, '%.2f  %.2f' % (analysis_series.mean(),
                                                  analysis_series.std()))


def abundances_to_df(pair, abundance_1, abundance_2):
    df = pd.DataFrame.from_dict({pair[:4]: abundance_1,
                                pair[4:]: abundance_2})

    return df


def get_characteristic_abundance_df(save_loc):
    out_df = pd.DataFrame()
    stats_df = pd.DataFrame()
    for pair in ['hisDgltA', 'hisDgltB', 'hisDpyrC']:
        strain_1_name = pair[:4]
        strain_2_name = pair[4:]

        def _drop_column(name):
            if 'AUX' not in name:
                return True
            if 'A0' in name or 'A1 ' in name or 'A7' in name:
                return True
            if 'I40' in name or 'I30' in name:
                return True
            return False

        df_raw = \
            pd.read_csv('%s/resequencing_data/AUX %s %s Mutation Table.csv' %
                        (resource_dir, strain_1_name, strain_2_name))

        df_raw.set_index('Position', inplace=True)
        df = df_raw[[i for i in df_raw.columns if not _drop_column(i)]]

        df = _fill_missing_df_entries(df, pair)

        def _get_average_mutation_fraction(strain_name):
            """Some characteristic mutations are more difficult for breseq
            to pick up so they appear as being present at 0%. In these cases
            defer to the value of the other characteristic mutation"""
            strain_fraction = np.zeros(len(df.columns))
            num_nonzero_array = np.zeros(len(df.columns))
            for characteristic in strain_to_characteristic[strain_name]:
                strain_fraction += df.loc[characteristic].values
                num_nonzero_array += (1. * df.loc[characteristic].values > 0)

            return strain_fraction / num_nonzero_array

        strain_1_abundances = _get_average_mutation_fraction(strain_1_name)
        strain_2_abundances = _get_average_mutation_fraction(strain_2_name)

        abundance_df = abundances_to_df(pair, strain_1_abundances,
                                        strain_2_abundances).T
        abundance_df.columns = df.columns
        abundance_df = _add_ale_numbers_to_df(abundance_df)

        abundance_df.to_csv('%s/abundance_by_characteristic_%s.csv' %
                            (save_loc, pair))

        _get_abundance_statistics(abundance_df, strain_2_name, stats_df)

        out_df = \
            _append_rows_to_abundance_df(out_df, strain_1_name, strain_2_name,
                                         strain_1_abundances,
                                         strain_2_abundances)
    print(out_df)
    out_df.to_excel('%s/characteristic_abundances_strain1.xlsx' % save_loc)
    stats_df.to_csv('%s/characteristic_hisD_abundance_stats.csv' % save_loc)


def get_coverage_abundance_df(save_loc, alignment_loc):
    out_df = pd.DataFrame()
    stats_df = pd.DataFrame()

    seqs = SeqIO.read('%s/resequencing_data/CP009273_1.gb' % resource_dir,
                      'gb')

    a = {}
    for feat in seqs.features:
        if feat.type != 'CDS':
            continue
        a[feat.qualifiers['gene'][0]] = \
            '%i-%i' % (feat.location.start, feat.location.end)

    i = 0
    for pair in ['hisD_gltA', 'hisD_gltB', 'hisD_pyrC']:
        df = pd.DataFrame()
        strain_1_name = pair.split('_')[0]
        strain_2_name = pair.split('_')[1]

        for file in glob.glob(alignment_loc +
                              'aux_%s/breseq/ale/*/data/reference.bam' % pair):
            ale_name = file.split('ale/')[1].split('/data')[0]
            summary_df = \
                pd.read_html(alignment_loc +
                             'aux_%s/breseq/ale/%s/output/summary.html'
                             % (pair, ale_name))
            fit_mean = float(summary_df[2][4][1])  # position of fit mean value

            # Filter starting strain and the strains w/ contamination
            if ale_name.split('-')[0] in ['1', '7', '0']:
                continue
            # Filter out endpoint strains
            if ale_name.split('-')[2] in ['30', '40']:
                continue

            cov_dict = {}
            for gene in ['hisD', 'pyrC', 'gltA', 'prpC', 'gdhA', 'gltB']:
                samfile = pysam.AlignmentFile(file, 'rb')
                read_cov = []
                start, stop = a[gene].split('-')
                for pileupcolumn in samfile.pileup('CP009273', int(start) + 30,
                                                   int(stop) - 30):
                    read_cov.append(pileupcolumn.n)
                samfile.close()
                cov_dict[gene] = np.array(read_cov).mean()

            ko1, ko2 = pair_to_kos[pair]
            avg_cov_ko1 = avg_cov_ko2 = 0
            i1 = i2 = 0
            for ko in ko1:
                print(ko, cov_dict[ko])
                i1 += 1
                avg_cov_ko1 += cov_dict[ko]
            print('--------------------------------')
            for ko in ko2:
                print(ko, cov_dict[ko])
                i2 += 1
                avg_cov_ko2 += cov_dict[ko]
            if np.isnan(avg_cov_ko1) and np.isnan(avg_cov_ko2):
                continue

            ale_info = ale_name.split('-')[:3]
            # abundance of strain is based on the percent of fit mean coverage
            # occupied by partner strain's KOed genes
            df.loc[i, strain_1_name] = avg_cov_ko2 / i2 / fit_mean
            df.loc[i, strain_2_name] = avg_cov_ko1 / i1 / fit_mean
            df.loc[i, 'Ale'] = float(ale_info[0])
            df.loc[i, 'Flask'] = float(ale_info[1])
            df.loc[i, 'Isolate'] = float(ale_info[2])

            i += 1

        abundance_df = \
            df.set_index(['Ale', 'Flask', 'Isolate']).T.sort_index(axis=1)
        abundance_df.to_csv('%s/abundance_by_coverage_%s.csv' %
                            (save_loc, pair.replace('_', '')))

        _get_abundance_statistics(abundance_df, strain_2_name, stats_df)

        out_df = \
            _append_rows_to_abundance_df(out_df, strain_1_name, strain_2_name,
                                         abundance_df.loc[strain_1_name].values,
                                         abundance_df.loc[strain_2_name].values)

    out_df.to_excel('%s/coverage_abundances_strain1.xlsx' % save_loc)
    stats_df.to_csv('%s/coverage_hisD_abundance_stats.csv' % save_loc)
