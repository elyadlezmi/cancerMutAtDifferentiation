import pandas as pd
import os
import projectFuncs as pf


VAR_TABLES = "/media/elyad/HDD/NeuralDiff/varTables"
STAR_COUNTS = "/media/elyad/HDD/NeuralDiff/counts"
DATA_TABLE = './CuratedMetaData.csv'
COSMIC = '/media/elyad/HDD/NeuralDiff/CosmicMutantExportCensus.tsv.gz'
BED_GRAPHS = "/media/elyad/HDD/NeuralDiff/bedGraphs"

MIN_COVERAGE = 10
MIN_ALT = 5
MIN_AF = 0
MIN_SNPS = 0
MIN_TUMORS = 10
MIN_SAMPLES_PER_LINE = 2
MIN_ALIGNED_READS = 2000000


# read metadata
data = pd.read_csv(DATA_TABLE).drop_duplicates(subset='Accession').fillna('NA')
data = data.loc[data['Accession'].isin([x.rstrip('Log.final.out') for x in os.listdir(STAR_COUNTS)])].reset_index(drop=True)
data = data.loc[data['Accession'].isin([x.rstrip('SNP.coverage') for x in os.listdir(BED_GRAPHS)])].reset_index(drop=True)
print(f"Remaining number of samples: {len(data)}")

# add stats
data = pf.meta_data_stats(data, vartable_path=VAR_TABLES, star_path=STAR_COUNTS,
                          min_coverage=MIN_COVERAGE, min_aligned_reads=MIN_ALIGNED_READS, min_alt=MIN_ALT)
print(f"Remaining number of samples: {len(data)}")

# analyze genotype, assign cell line ID
data = pf.annotate_cell_lines(data, vartable_path=VAR_TABLES, bed_path=BED_GRAPHS, min_coverage=15)
data.to_csv('metaDataWithCellLineID.csv')
data = pd.read_csv('./metaDataWithCellLineID.csv')

# filter samples with low SNP coverage, bad ratio or cell lines that have only one sample
data, excluded = pf.filterMetaData(data, min_snps=MIN_SNPS, min_samples_per_line=MIN_SAMPLES_PER_LINE)
print(f"final number of samples: {len(data)}, excluded samples: {len(excluded)}")

# analyze for cancer mutations and aggregate all muations
cosmic = pd.read_csv(COSMIC, sep='\t', compression='gzip', encoding='latin1')

data, mutations = pf.aggregate_mutations(data, cosmic, min_coverage=MIN_COVERAGE, min_tumors=MIN_TUMORS, min_af=MIN_AF, vartable_path=VAR_TABLES, min_alt=MIN_ALT)
excluded, excluded_mutations = pf.aggregate_mutations(excluded, cosmic, min_coverage=MIN_COVERAGE, min_tumors=MIN_TUMORS, min_af=MIN_AF, vartable_path=VAR_TABLES, min_alt=MIN_ALT)

# add info whether mutations are somatic mutations (acquired in culture)
mutations = pf.is_somatic(mutations, data, min_coverage=MIN_COVERAGE, path_to_bed=BED_GRAPHS)

# write excluded samples mutations
excluded.to_csv('finalExcludedMetaData.csv')
with open('excludedMutations.csv', 'w') as f:
    f.write(f"Samples: {len(excluded['Accession'])} Studies: {len(set(excluded['Study']))} Cell lines:\
            {len(set(excluded['CellLineID']))} params: {MIN_COVERAGE} {MIN_ALT} {MIN_AF} {MIN_SNPS} \
            {MIN_TUMORS} {MIN_SAMPLES_PER_LINE}\n")
    excluded_mutations.to_csv(f, index=False)

# write tables 1.All mutations 2. Agg by COSMIC_ID 3. Agg by gene
data.to_csv('finalMetaData.csv')
with open('Mutations.csv', 'w') as f:
    f.write(f"Samples: {len(data['Accession'])} Studies: {len(set(data['Study']))} Cell lines:\
            {len(set(data['CellLineID']))} params: {MIN_COVERAGE} {MIN_ALT} {MIN_AF} {MIN_SNPS}\
            {MIN_TUMORS} {MIN_SAMPLES_PER_LINE}\n")
    mutations.to_csv(f, index=False)

summary = pf.groupByMutation(mutations)
summary.to_csv('MutationSummary.csv')
pf.groupByGene(summary).to_csv('MutationSummarizedPerGene.csv')

# acquired mutations
acquired = pf.acquired_mutations(mutations, data, min_coverage=MIN_COVERAGE, bed_path=BED_GRAPHS)

# write tables > 1.All mutations 2. Agg by COSMIC_ID 3. Agg by gene
acquired.to_csv('AcquiredMutations.csv', index=False)
summary = pf.groupByMutation(acquired)
summary.to_csv('AcquiredMutationSummary.csv')
pf.groupByGene(summary).to_csv('AcquiredMutationSummarizedPerGene.csv')






