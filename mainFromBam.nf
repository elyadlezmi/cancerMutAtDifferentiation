#!/usr/bin/env nextflow

// pipeline for identification of cancer related mutations

prefix = params.prefix


process markDuplicates {
    
    memory '32GB'
    publishDir "$launchDir", mode: 'copy'
	
    input:
    path bam from Channel.fromPath("${prefix}Filtered.bam") 
    
    output:
    path "${prefix}.bam" into markedDuplicates
    path "${prefix}.bam" into markedDuplicatesForIndex
        
    """
    /gatk-4.1.3.0/gatk MarkDuplicates --I $bam --O ${prefix}.bam --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT --M ${prefix}marked_dup_metrics.txt
    """ 
}

process indexBam {
    
    memory '16GB'
    publishDir "$launchDir", mode: 'copy'
	
    input:
    path bam from markedDuplicatesForIndex
    
    output:
    path "${prefix}.bam.bai" into indexedForSplitN
    """
    samtools index $bam
    """ 
}
   
   
process splitNcigar {
    
    memory '16GB'
    
    input:
    path bam from markedDuplicates
    path index from indexedForSplitN
    path intervals from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz")
    path intIndex from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz.tbi")
    path reference_genome from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa")
    path index from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa.fai")
    path dict from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.dict")
    
    output:
    path "${prefix}splitN.bam" into splitN 

    """
    /gatk-4.1.3.0/gatk SplitNCigarReads -L $intervals -R $reference_genome -I $bam -O ${prefix}splitN.bam
    """ 
}
    
    
process addGroups {
    
    memory '16GB'
    
    input:
    path bam from splitN
    
    output:
    path "${prefix}.grouped.bam" into grouped4BQSR
    path "${prefix}.grouped.bai" into grouped4BQSRindex
    path "${prefix}.grouped.bam" into grouped4applyBQSR
    path "${prefix}.grouped.bai" into grouped4applyBQSRindex
    """
    /gatk-4.1.3.0/gatk AddOrReplaceReadGroups --CREATE_INDEX true --I $bam --O ${prefix}.grouped.bam --RGID rnasq --RGLB lb --RGPL illumina --RGPU pu --RGSM $prefix
    """  
}

    
process baseQualityRecalibration {
    
    memory '16GB'
    
    input:
    path bam from grouped4BQSR
    path index from grouped4BQSRindex
    path intervals from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz")
    path intIndex from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz.tbi")
    path reference_genome from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa")
    path index from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa.fai")
    path dict from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.dict")
    path dbSNP from Channel.fromPath("$projectDir/data/dbSNPbuild154Renamed.vcf.gz")
    path dbSNPindex from Channel.fromPath("$projectDir/data/dbSNPbuild154Renamed.vcf.gz.tbi")
    
    output:
    path "${prefix}.recal_data.table" into recalTable 

    """
    /gatk-4.1.3.0/gatk BaseRecalibrator -L $intervals -I $bam --use-original-qualities --disable-sequence-dictionary-validation true -R $reference_genome --known-sites $dbSNP -O ${prefix}.recal_data.table
    """ 
}
    
    
process applyBQSR {
    
    memory '16GB'
    
    input:
    path table from recalTable
    path bam from grouped4applyBQSR
    path bai from grouped4applyBQSRindex
    path intervals from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz")
    path intIndex from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz.tbi")
    path reference_genome from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa")
    path index from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa.fai")
    path dict from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.dict")
    
    output:
    path "${prefix}.recal_output.bam" into recalibrated
    path "${prefix}.recal_output.bai" into recalibratedIndex  

    """
    /gatk-4.1.3.0/gatk ApplyBQSR -L $intervals -R $reference_genome -I $bam --use-original-qualities --add-output-sam-program-record --bqsr-recal-file $table -O ${prefix}.recal_output.bam
    """ 
}



process callVariants {
    
    if (params.keepInter == true) {
        publishDir "$launchDir", mode: 'copy'}
    memory '16GB'
    
    input:
    path bam from recalibrated
    path bai from recalibratedIndex
    path intervals from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz")
    path intIndex from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz.tbi")
    path reference_genome from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa")
    path index from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa.fai")
    path dict from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.dict")
    
    output:
    file "${prefix}.output.vcf.gz" into variants
    file "${prefix}.output.vcf.gz.tbi" into variantsIndex 

    """
    /gatk-4.1.3.0/gatk HaplotypeCaller -L $intervals -R $reference_genome -I $bam -O ${prefix}.output.vcf.gz --dont-use-soft-clipped-bases --pcr-indel-model AGGRESSIVE
    """ 
}

    
process hardFilter {
    
    memory '16GB'
    
    input:
    path vcf from variants
    path vcfIndex from variantsIndex
    path reference_genome from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa")
    path index from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa.fai")
    path dict from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.dict")
    
    output:
    path "${prefix}.hardfilter.vcf.gz" into hardfilterd


    """
    /gatk-4.1.3.0/gatk VariantFiltration --R $reference_genome --V $vcf --window 35 --cluster 3 --filter-name "FS" --filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0" -O ${prefix}.hardfilter.vcf.gz
    """ 
}    

    
process filterVariants {
    
    cpus params.cpu
    
    input:
    path hardfilterd
   
    output:
    path "${prefix}filtered.vcf.gz" into filterd

    """
    bcftools view --threads $params.cpu -i 'FILTER="PASS"' --output-type z --output-file ${prefix}filtered.vcf.gz ${prefix}.hardfilter.vcf.gz
    """ } 

process IndexfilteredVariants {
        
    input:
    path filterd
   
    output:
    path "${prefix}filtered.vcf.gz" into filterd4annotation
    path "${prefix}filtered.vcf.gz.tbi" into filteredIndex  

    """
    tabix "${prefix}filtered.vcf.gz"
    """ } 
      
process geneAnnotation {
    
    cpus params.cpu
    
    input:
    path filterd4annotation
    path intervals from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz")
    path header from Channel.fromPath("$projectDir/data/header.txt")
    path filteredIndex
   
    output:
    path "${prefix}named.vcf.gz" into gene
    path "${prefix}named.vcf.gz.tbi" into geneIndex  

    """
    bcftools annotate --threads $params.cpu -a $intervals -h $header -c CHROM,FROM,TO,Gene --output-type z --output ${prefix}named.vcf.gz ${prefix}filtered.vcf.gz
    tabix ${prefix}named.vcf.gz
    """ } 

process snpAnnotation {
    
    cpus params.cpu
    
    input:
    path gene
    path geneIndex
    path dbSNP from Channel.fromPath("$projectDir/data/dbSNPbuild154Renamed.vcf.gz")
    path dbSNPindex from Channel.fromPath("$projectDir/data/dbSNPbuild154Renamed.vcf.gz.tbi")
    
    output:
    path "${prefix}dbSNP.vcf.gz" into snp
    path "${prefix}dbSNP.vcf.gz.tbi" into snpIndex  

    """
    bcftools annotate --threads $params.cpu -a $dbSNP -c INFO/RS,INFO/COMMON --output-type z --output ${prefix}dbSNP.vcf.gz ${prefix}named.vcf.gz
    tabix ${prefix}dbSNP.vcf.gz
    """ } 

process cosmicAnnotation {
   
    cpus params.cpu
       
    input:
    path snp
    path snpIndex
    path cosmic_vcf from Channel.fromPath("$projectDir/data/CosmicCodingMutsRenamed.vcf.gz")
    path cosmicIndex from Channel.fromPath("$projectDir/data/CosmicCodingMutsRenamed.vcf.gz.tbi")
    
    output:
    path "${prefix}.annotated.vcf.gz" into cosmic
    path "${prefix}.annotated.vcf.gz.tbi" into cosmicIndex  

    """
    bcftools annotate --threads $params.cpu -a $cosmic_vcf -c ID,INFO/CNT --output-type z --output ${prefix}.annotated.vcf.gz ${prefix}dbSNP.vcf.gz
    tabix ${prefix}.annotated.vcf.gz
    """
    } 
    
process variantTable {
    
    publishDir "$launchDir", mode: 'copy'
    
    input:
    path cosmic
    path cosmicIndex
   
    output:
    path "${prefix}varTable.tsv" into table  

    """
    bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/Gene\t%ID\t%INFO/CNT\t%INFO/RS\t%INFO/COMMON\t[%AD]\n' ${prefix}.annotated.vcf.gz > ${prefix}varTable.tsv
    """
    } 
