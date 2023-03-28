#!/usr/bin/python3
#SBATCH --mem=4G
#SBATCH -c1
#SBATCH --time=48:00:00

import pandas as pd
import os
import time

PROJECT_DIR = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations"
DATA_TABLE = '/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/CuratedMetaData.csv'
BAMS = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/undupBams"
FASTQS = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/fastqs"
COUNTS = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/counts"
VAR_TABLES = "/sci/labs/nissimb/elyad/icore-data/NeuralDiffMutations/varTables"
STAR_INDEX="/sci/labs/nissimb/elyad/icore-data/STARindex"

data = pd.read_csv(DATA_TABLE)  # df of accession only
data = data.loc[data['Accession'].isin([x.rstrip('varTable.tsv') for x in os.listdir(VAR_TABLES)])]

#time.sleep(3600)
n = 0
done = []
for i, row in data.iloc[::].iterrows():
	
	sample = row['Accession']
	
	# done
	if os.path.isfile(f"{COUNTS}/{sample}ReadsPerGene.out.tab"):
		n += 1
		done.append(sample)
		continue
	
	# recently done
	if os.path.isfile(f"{PROJECT_DIR}/workDir/{sample}/STARdone.txt"):
		os.system(f"mv {PROJECT_DIR}/workDir/{sample}/{sample}ReadsPerGene.out.tab {COUNTS}")
		os.system(f"mv {PROJECT_DIR}/workDir/{sample}/{sample}Log.final.out {COUNTS}")
		os.system(f"rm -rf {PROJECT_DIR}/workDir/{sample}")
		continue
	
	if not os.path.isfile(f"{PROJECT_DIR}/workDir/{sample}/done.txt"):
		continue
	
	os.chdir(f"{PROJECT_DIR}/workDir/{sample}")
	os.system("rm -rf done.txt")
	
	with open(f"{sample}count.sh", 'w') as f:
		f.write(f"#!/usr/bin/bash\n")
		f.write("#SBATCH--mem=32G\n")
		f.write("#SBATCH-c8\n")
		f.write("#SBATCH--time=24:00:00\n")
		f.write(f"module load STAR/2.7.10a\n""")
		if row['LibraryLayout'] == 'single' or row['LibraryLayout'] == 'SINGLE':
			f.write(f"STAR --runThreadN 8 --genomeDir {STAR_INDEX} --readFilesIn {sample}.fastq.gz --outFileNamePrefix {sample} --outSAMtype BAM Unsorted --readFilesCommand zcat --quantMode GeneCounts\n")
		else:
			f.write(f"STAR --runThreadN 8 --genomeDir {STAR_INDEX} --readFilesIn {sample}_1.fastq.gz {sample}_2.fastq.gz --outFileNamePrefix {sample} --outSAMtype BAM Unsorted --readFilesCommand zcat --quantMode GeneCounts\n")

		f.write(f"if ( test -f {sample}ReadsPerGene.out.tab ) && ( test -f {sample}Log.final.out ); then echo XXXXXXX > STARdone.txt; fi\n")
		
	os.system(f"sbatch {sample}count.sh")
	
print(n)
print(len(data))
data = data.loc[~data['Accession'].isin(done)]
data.to_csv(f'{PROJECT_DIR}/data2align.csv')

		
		

