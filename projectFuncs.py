import os
import re
import numpy as np
import pandas as pd
from natsort import natsorted, natsort_keygen

def cleanVarTable(varTablePath):

    df = pd.read_csv(varTablePath, sep='\t')
    df.columns = ['CHROM', 'POS', 'REF', 'ALT', 'Gene', 'COSMIC_ID', 'COSMIC_CNT', 'RS(dbSNP)', 'COMMON', 'AD']
    df[['Ref #', 'Alt #']] = df['AD'].str.split(",", expand=True).iloc[:, [0, 1]]
    df = df.astype({'Ref #': int, 'Alt #': int})
    df['Reads'] = df[['Ref #', 'Alt #']].sum(axis=1)
    df['AF'] = df['Alt #'].divide(df['Reads'], fill_value=1)
    del df['AD']

    return df


def varToMut(df, cosmic, min_reads=10, min_tumors=10, min_af=0, min_alt=2):

    cosmic = cosmic[['GENOMIC_MUTATION_ID', 'Mutation AA', 'Mutation Description', 'FATHMM score', 'Tier', 'FATHMM prediction']]

    df = df.loc[df['Reads'] >= min_reads]

    df = df.loc[df['COSMIC_CNT'] != '.']
    df = df.astype({'COSMIC_CNT': int})
    df = df.loc[(df['COSMIC_CNT'] >= min_tumors) & (df['COMMON'] != '1') & (df['AF'] > min_af) & (df['Alt #'] >= min_alt)]

    df = df.merge(cosmic, left_on='COSMIC_ID', right_on='GENOMIC_MUTATION_ID')

    df = df.loc[(df['Tier'] == 1) & (~df['Mutation Description'].isin(['Substitution - coding silent', 'Unknown'])) & (
            df['FATHMM prediction'] == 'PATHOGENIC')]
    df.drop_duplicates(subset='COSMIC_ID', inplace=True)

    del df['Tier']
    del df['FATHMM prediction']
    del df['GENOMIC_MUTATION_ID']

    return df


def getPositions(df):
    df = df.iloc[:, [0,1]].astype(str)
    df = df.iloc[:, 0] + '_' + df.iloc[:, 1]
    return list(df)


def get_covered_pos(sample, bed_path, min_coverage):
    # if os.path.getsize(f"{bed_path}/{sample}SNP.coverage") != 0:
    pos = pd.read_csv(f"{bed_path}/{sample}SNP.coverage", sep='\t', header=None)
    return getPositions(pos.loc[pos.iloc[:, 2] >= min_coverage])


def get_common_snps(sample, vartable_path, min_coverage):
    df = cleanVarTable(f"{vartable_path}/{sample}varTable.tsv")
    df = df.loc[(df['Reads'] >= min_coverage) & (df['COMMON'] == '1') & (df['CHROM'] != 'chrM')].reset_index(drop=True)
    df['position'] = getPositions(df)
    return df


def get_covered_snps(df, covered_pos):
    df = df.loc[df['position'].isin(covered_pos)]
    df = df.iloc[:, [0, 1, 3]].astype(str)
    df = df.iloc[:, 0] + '_' + df.iloc[:, 1] + '_' + df.iloc[:, 2]
    return list(df)


def is_same_line(snp1, snp2):
    snp1, snp2 = set(snp1), set(snp2)
    a, b, c = len(snp1.intersection(snp2)), len(snp1.difference(snp2)), len(snp2.difference(snp1))
    if a+b+c == 0:
        return False
    else:
        return (b+c)/(a+b+c) < 0.3


def annotate_cell_lines(data, vartable_path, bed_path, min_coverage):

    # sort to start from best samples
    data = data.sort_values(by="n_snps", ascending=True)

    samples_left = [
        (sample, get_covered_pos(sample, bed_path, min_coverage), get_common_snps(sample, vartable_path, min_coverage))
        for sample in data['Accession']]

    # rearrange for faster run
    cell_lines, n = {}, 0

    while len(samples_left) > 0:

        print(f"samples left: {len(samples_left)}")
        n += 1
        f_sample, f_coverage, f_snps = samples_left.pop()
        cell_lines[f"CL{n}"] = [f_sample]
        samples_to_remove = []
        for i, value in enumerate(samples_left):
            sample, coverage, snps = value
            shared = natsorted(list(set(coverage).intersection(f_coverage)))
            g1, g2 = get_covered_snps(f_snps, shared), get_covered_snps(snps, shared)
            if is_same_line(g1, g2):
                cell_lines[f"CL{n}"].append(sample)
                samples_to_remove.append(i)

            # print(data.loc[data['Accession'].isin([sample, f_sample]), ['Accession', 'CellLine', 'Study', 'CellState']])
        for i in reversed(samples_to_remove):
            samples_left.pop(i)

    s, l = [], []
    for key in cell_lines.keys():
        for val in cell_lines[key]:
            s.append(val)
            l.append(key)

    data = pd.merge(data, pd.DataFrame({'Accession': s, 'CellLineID': l}), on='Accession')
    return data


def aggregate_mutations(data, cosmic, min_coverage, min_tumors, min_af, min_alt, vartable_path):
    df, n_mut, muts, af, n_Tpoints = pd.DataFrame(), [], [], [], []
    temp_df = data.groupby('Study')['TimePoint'].nunique()
    for i, row in data.iterrows():
        p = cleanVarTable(f"{vartable_path}/{row['Accession']}varTable.tsv")
        p = varToMut(p, cosmic, min_reads=min_coverage, min_tumors=min_tumors, min_af=min_af, min_alt=min_alt)
        n_mut.append(len(p))
        muts.append(' '.join(list(p['Gene'])))
        af.append(' '.join(list(p['AF'].astype(str))))
        p[['Accession', 'Study', 'CellState', 'CellLine', 'SampleType', 'TimePoint', 'iPSES', 'CellLineID']] = \
            row[['Accession', 'Study', 'CellState', 'CellLine',  'SampleType', 'TimePoint', 'iPSES', 'CellLineID']]
        df = pd.concat([df, p])
        n_Tpoints.append(temp_df[row['Study']])
    data['n_mut'] = n_mut
    data['muts'] = muts
    data['AF'] = af
    data['n_TimePoints'] = n_Tpoints
    return data, df


def isPosExpressedInSamples(pos, samples, threshold, path_to_bed):
    is_expressed = True
    for sample in samples:
        coverage = pd.read_csv(f"{path_to_bed}/{sample}Var.coverage", sep='\t', header=None)
        condition = coverage.iloc[:, 2] >= threshold
        coverage = coverage.loc[condition]
        is_expressed = pos in getPositions(coverage)
        if not is_expressed:
            break

    return is_expressed


def groupByMutation(df):
    df = df.groupby(by='COSMIC_ID').agg({'Gene': 'first', 'CHROM': 'first', 'POS': 'first', 'REF': 'first', 'ALT': 'first',
                                         'COSMIC_CNT': 'first', 'RS(dbSNP)': 'first',  'COMMON': 'first', 'Accession': list,
                                         'Study': set, 'CellLineID': set, 'TimePoint': list})
    df[['N_samples', 'N_studies', 'N_lines']] = df[['Accession', 'Study', 'CellLineID']].applymap(len)

    return df.sort_values(by='N_studies', ascending=False)


def groupByGene(mutationSummary):
    mutationSummary['N_mutations'] = mutationSummary['Gene']
    df = mutationSummary.groupby(by='Gene').agg(
        {'Gene': 'first', 'Accession': list, 'Study': list, 'CellLineID': list, 'N_mutations': 'count'})
    df[['Accession', 'Study', 'CellLineID']] = df[['Accession', 'Study', 'CellLineID']].applymap(lambda x: set([i for j in x for i in j]))
    df[['N_samples', 'N_studies', 'N_lines']] = df[['Accession', 'Study', 'CellLineID']].applymap(len)

    return df.sort_values(by='CellLineID', ascending=False)


def acquired_mutations(mutations, data, min_coverage, bed_path):
    acquired = []
    for i, row in mutations.loc[mutations['CellState'] == 'differentiated'].iterrows():
        study, cell_line, cosmic_id = str(row['Study']), str(row['CellLineID']), str(row['COSMIC_ID'])
        t_data = data.loc[(data['Study'] == study) & (data['CellLineID'] == cell_line)]
        early_time = natsorted(list(t_data['TimePoint']))[0]

        if row['TimePoint'] == early_time:
            continue

        early_samples = t_data.loc[t_data['TimePoint'] == early_time, 'Accession']

        # if len(early_samples) == 0:
        #     continue

        if cosmic_id in list(mutations.loc[mutations['Accession'].isin(early_samples), 'COSMIC_ID']):
            continue

        # is expressed in hESCs
        pos = str(row['CHROM']) + '_' + str(row['POS'])
        if isPosExpressedInSamples(pos, early_samples, min_coverage, bed_path):
            acquired.append(list(row))

    return pd.DataFrame(acquired, columns=mutations.columns)


def is_somatic(mutations, data, path_to_bed, min_coverage):

    mutations['mut_line'] = mutations['COSMIC_ID'] + '_' + mutations['CellLineID']
    m = mutations.drop_duplicates(subset='mut_line')
    mut_line, isit_somatic = [], []

    for i, row in m.iterrows():

        pos = str(row['CHROM']) + '_' + str(row['POS'])
        samples = list(data.loc[(data['CellLineID'] == row['CellLineID']) & (~data['Accession'].isin(mutations['Accession'])), 'Accession'])

        somatic = False
        for sample in samples:
            coverage = pd.read_csv(f"{path_to_bed}/{sample}Var.coverage", sep='\t', header=None)
            condition = coverage.iloc[:, 2] >= min_coverage
            coverage = coverage.loc[condition]

            if pos in getPositions(coverage):
                somatic = True
                break

        mut_line.append(row['mut_line'])
        isit_somatic.append(somatic)

    res = pd.DataFrame({'mut_line': mut_line, 'is_somatic': isit_somatic})
    mutations = pd.merge(mutations, res, on='mut_line')
    del mutations['mut_line']

    return mutations


def meta_data_stats(data, vartable_path, star_path, min_coverage, min_aligned_reads, min_alt):
    n_reads, lAR, n_vars, n_snps = [], [], [], []
    for sample in data['Accession']:
        df = cleanVarTable(f"{vartable_path}/{sample}varTable.tsv")
        df = df.loc[df['Reads'] >= min_coverage]
        df = df.loc[df['Alt #'] >= min_alt]

        n_vars.append(len(df))
        n_snps.append(len(df.loc[df['RS(dbSNP)'] != '.']))

        if os.path.isfile(f"{star_path}/{sample}Log.final.out"):
            f = open(f"{star_path}/{sample}Log.final.out")
            t = f.readlines()
            p = re.findall(r"\d+\.\d+", t[9])
            n = re.findall(r"\d+", t[5])
            lAR.append(float(p[0]))
            n_reads.append(int(n[0]))

    data = pd.concat([data.reset_index(drop=True),
                      pd.DataFrame({'n_reads': n_reads, 'AlignmentRate': lAR, 'n_vars': n_vars, 'n_snps': n_snps})], axis=1)
    data = data.loc[data['n_reads']*data['AlignmentRate']/100 >= min_aligned_reads]

    return data.reset_index(drop=True)


def filterMetaData(data, min_snps, min_samples_per_line):
    before = data.loc[:]
    data = data.loc[data['n_snps'] >= min_snps]


    for cl in set(data['CellLineID']):
        if len(data.loc[data['CellLineID'] == cl, 'Accession']) < min_samples_per_line:
            data = data.loc[data['CellLineID'] != cl]

    excluded = set(before['Accession']).difference(set(data['Accession']))
    excluded = before.loc[before['Accession'].isin(excluded)]

    return data.reset_index(drop=True), excluded.reset_index(drop=True)
