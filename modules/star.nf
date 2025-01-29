//Indexing the genome with STAR
process STAR_INDEX {
    conda "envs/map.yaml"
    label 'high'
    publishDir "$params.outdir/STAR/", mode:'copy'
    
    output:
    path("Indicies")
    
    when:
    !params.isLong
    
    script:
    """
    mkdir -p Indicies
    STAR --runThreadN 14 --runMode genomeGenerate --genomeDir Indicies --genomeFastaFiles ${params.genome} --sjdbGTFfile ${params.annotation} --sjdbOverhang ${params.overhang} --limitGenomeGenerateRAM 60000000000
    """
}

//Mapping filtered short single/paired-end reads to genome generating a coordinate-sorted bam file
process STAR_MAP {
    conda "envs/map.yaml"
    label 'high'
    tag "on $read_id"
    publishDir "$params.outdir/STAR/", mode:'copy'
    
    input:
    path indices
    tuple val(read_id), path(reads)
    
    when:
    !params.isLong
    
    output:
    path("map_${read_id}")
    tuple val(read_id), path("map_${read_id}/${read_id}_Aligned.sortedByCoord.out.bam"), emit: sorted_bam
    path("map_${read_id}/${read_id}_Log.final.out"), emit: log_final
    tuple val(read_id), path("map_${read_id}/${read_id}.bam.bai"), emit: index
    
    
    script:
    """
    mkdir -p map_${read_id}
    if [[ ${params.isSingle} == true ]]; then
        STAR --runThreadN 28 --genomeDir ${indices} --readFilesIn ${reads}  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix map_${read_id}/${read_id}_
    else
        STAR --runThreadN 28 --genomeDir ${indices} --readFilesIn ${reads[0]} ${reads[1]}  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix map_${read_id}/${read_id}_
    fi
    samtools index map_${read_id}/${read_id}_Aligned.sortedByCoord.out.bam map_${read_id}/${read_id}.bam.bai
    """
}
