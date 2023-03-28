#!/usr/bin/env nextflow

// pipeline for identification of cancer related mutations

prefix = params.prefix

if ( params.fastq2 == false ) {
    
    process trimmomaticSE {
        
        cpus params.cpu
	    memory '16GB'
        
        input:
        path fastq from Channel.fromPath("$params.fastq") 
        path adapters from Channel.fromPath("$projectDir/data/CommonAdapters.fa")
        
        output: 
        path "${prefix}.trimmed.fastq.gz" into trimmed, trimmed4mouse
        
        """
        java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $params.cpu $fastq ${prefix}.trimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:26
        """ 
    }
      
             
    process alignToHumanGenomeSE {
          
        cpus params.cpu
		memory '64GB'
        
        input:
        path fastq from trimmed
        path genomeDir from Channel.fromPath("$projectDir/data/GRCh", type: 'dir')
        
        output:
        path "${prefix}Aligned.sortedByCoord.out.bam" into aligned2human
        val 'foo' into bar
        
        """
        STAR --runThreadN $params.cpu --genomeDir $genomeDir --readFilesIn $fastq --outFileNamePrefix ${prefix} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMattributes NM --twopassMode Basic --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1
        """
    }


    process alignToMouseGenomeSE {
        
        cpus params.cpu
		memory '64GB'
           
        input:
        path fastq from trimmed4mouse
        path genomeDir from Channel.fromPath("$projectDir/data/GRCm", type: 'dir')
        val 'foo' from bar
        
        output:
        path "${prefix}GRCmAligned.sortedByCoord.out.bam" into aligned2mouse
        
        when:
        params.filterMouse == true
        
        """
        STAR --runThreadN $params.cpu --genomeDir $genomeDir --readFilesIn $fastq --outFileNamePrefix ${prefix}GRCm --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMattributes NM --twopassMode Basic --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1
        """ 
    }

} else {

    process trimmomaticPE {
        
        cpus params.cpu
        memory '16GB'
	
        input:
        path fastq1 from Channel.fromPath("$params.fastq") 
        path fastq2 from Channel.fromPath("$params.fastq2") 
        path adapters from Channel.fromPath("$projectDir/data/CommonAdapters.fa")
         
        output: 
        path "${prefix}_1.trimmed.fastq.gz" into trimmed1
        path "${prefix}_2.trimmed.fastq.gz" into trimmed2
        path "${prefix}_1.trimmed.fastq.gz" into trimmed4mouse1
        path "${prefix}_2.trimmed.fastq.gz" into trimmed4mouse2
        
        """
        java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $params.cpu $fastq1 $fastq2 ${prefix}_1.trimmed.fastq.gz unpaired1.fastq.gz ${prefix}_2.trimmed.fastq.gz unpaired2.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:26
        """ 
    }


    process alignToHumanGenomePE {
        
        cpus params.cpu
		memory '64GB'
           
        input:
        path fastq1 from trimmed1
        path fastq2 from trimmed2
        path genomeDir from Channel.fromPath("$projectDir/data/GRCh", type: 'dir')
        
        output:
        path "${prefix}Aligned.sortedByCoord.out.bam" into aligned2human
        
        """
        STAR --runThreadN $params.cpu --genomeDir $genomeDir --readFilesIn $fastq1 $fastq2 --outFileNamePrefix ${prefix} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMattributes NM --twopassMode Basic --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1
        """ 
    }
     
        
    process alignToMouseGenomePE {
        
        cpus params.cpu
		memory '64GB'
                
        input:
        path fastq1 from trimmed4mouse1
        path fastq2 from trimmed4mouse2
        path genomeDir from Channel.fromPath("$projectDir/data/GRCm", type: 'dir')
        
        output:
        path "${prefix}GRCmAligned.sortedByCoord.out.bam" into aligned2mouse
        
        when:
        params.filterMouse == true
        
        """
        STAR --runThreadN $params.cpu --genomeDir $genomeDir --readFilesIn $fastq1 $fastq2 --outFileNamePrefix ${prefix}GRCm --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMattributes NM --twopassMode Basic --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1
        """ 
        }
}

if ( params.filterMouse == true ) {
    
    process XenofilteR {
        
        publishDir "$launchDir", mode: 'copy'
        cpus params.cpu
		memory '64GB'
        
        input:
        path bam1 from aligned2human
        path bam2 from aligned2mouse
        
        output:
        file "Filtered_bams/${prefix}_Filtered.bam" into filteredBam
		
        """
        #!/usr/bin/Rscript --save
        library("XenofilteR")
        bp.param <- SnowParam(workers = $params.cpu, type = "SOCK")
        sample.list <- matrix(c('$bam1','$bam2'),ncol=2)
        output.names <- c('${prefix}')
        XenofilteR(sample.list, destination.folder = "./", MM_threshold = 8, bp.param = bp.param, output.names)
        """ 
    }
    
    
} else {

    process skipXenofilteR {
        
		if (params.keepInter == true) {
            publishDir "$launchDir", mode: 'copy'}
        
        input:
        path bam from aligned2human
        
        output:
        path bam into filteredBam 
        
        """
        echo "--- mouse read filtering was not performed ---"
        """
    }
}

