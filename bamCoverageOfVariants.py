#!/usr/bin/python3
# SBATCH --mem=32G
# SBATCH -c1
# SBATCH --time=48:00:00

import pandas as pd
import os


def cleanVarTable(varTablePath):
    """Arrange table"""
    df = pd.read_csv(varTablePath, sep='\t')
    df.columns = ['CHROM', 'POS', 'REF', 'ALT', 'Gene', 'COSMIC_ID', 'COSMIC_CNT', 'RS(dbSNP)', 'COMMON', 'AD']
    df[['Ref #', 'Alt #']] = df['AD'].str.split(",", expand=True).iloc[:, [0, 1]]
    df = df.astype({'Ref #': int, 'Alt #': int})
    df['Reads'] = df[['Ref #', 'Alt #']].sum(axis=1)
    df['AF'] = df['Alt #'].divide(df['Reads'], fill_value=1)
    del df['AD']

    return df


DATA_TABLE = '/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/CuratedMetaData.csv'
BED_GRAPPHS = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/bedGraphs"
VAR_TABLES = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/varTables"
BAMS = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/undupBams"
DBSNP = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/dbSNP154commonExome.bed"
ALL_VAR_POS = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/allVarPositions.bed"
GTF = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/gencode.v38.primary_assembly.annotation.gtf"
COUNTS = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/counts"

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

data = pd.read_csv(DATA_TABLE, usecols=['Accession', 'Study', 'CellState', 'CellLine'])  # df of accession only
data = data.loc[data['Accession'].isin([x.rstrip('varTable.tsv') for x in os.listdir(VAR_TABLES)])]

samples = data['Accession']

# create bed with all variant positions in all studies
#df, df_sorted = pd.DataFrame(columns=['CHROM', 'POS']), pd.DataFrame(columns=['CHROM', 'POS'])
#for sample in samples:
#    variants = cleanVarTable(f"{VAR_TABLES}/{sample}varTable.tsv")
#    df = pd.concat([df, variants[['CHROM', 'POS']]]).drop_duplicates()
#for chrom in chromosomes:
#    df_sorted = pd.concat([df_sorted, df.loc[df['CHROM'] == chrom].sort_values(by='POS')])
#df_sorted.to_csv('/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/allVarPositions.bed', header=None, index=False,
#                 sep='\t')

no_bam = []
for sample in samples:
    if not os.path.isfile(f"{BAMS}/{sample}.bam"):
        no_bam.append(sample)

with open('no_bams.txt', 'w') as f:
    for i in no_bam:
        f.write(i + '\n')
        samples.remove(sample)

for sample in samples:

    if not os.path.isfile(f"{BED_GRAPPHS}/{sample}Var.coverage"):
        with open(f"{sample}Bams1.txt", 'w') as f:
            f.write(f"{BAMS}/{sample}.bam\n")

        with open(f"{sample}Vars.sh", 'w') as f:
            f.write(f"#!/usr/bin/bash\n")
            f.write("#SBATCH--mem=8G\n")
            f.write("#SBATCH-c1\n")
            f.write("#SBATCH--time=24:00:00\n")
            f.write(
                f"samtools depth -b {ALL_VAR_POS} -f {sample}Bams1.txt | " + """awk -F"\\t" '$3 > 5  { print $1"\\t"$2"\\t"$3 }' > """ + f"{BED_GRAPPHS}/{sample}Var.coverage\n")
            f.write(f"rm -rf {sample}Vars.sh {sample}Bams1.txt")

        os.system(f"sbatch {sample}Vars.sh")

    if not os.path.isfile(f"{BED_GRAPPHS}/{sample}SNP.coverage"):
        with open(f"{sample}Bams2.txt", 'w') as f:
            f.write(f"{BAMS}/{sample}.bam\n")

        with open(f"{sample}SNPS.sh", 'w') as f:
            f.write(f"#!/usr/bin/bash\n")
            f.write("#SBATCH--mem=8G\n")
            f.write("#SBATCH-c1\n")
            f.write("#SBATCH--time=24:00:00\n")
            f.write(
                f"samtools depth -b {DBSNP} -f {sample}Bams2.txt | " + """awk -F"\\t" '$3 > 5  { print $1"\\t"$2"\\t"$3 }' > """ + f"{BED_GRAPPHS}/{sample}SNP.coverage\n""")
            f.write(f"rm -rf {sample}Bams2.txt {sample}SNPS.sh\n")

        os.system(f"sbatch {sample}SNPS.sh")
