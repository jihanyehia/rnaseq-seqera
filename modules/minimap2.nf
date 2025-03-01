//Building index for Minimap2
process MINIMAP2_INDEX {
    container 'community.wave.seqera.io/library/minimap2_picard_samtools_star_subread:a1a7ac39cfcdefee'
    memory '20 GB'
    publishDir "$params.outdir/Minimap2/Index", mode:'copy'
    
    input:
    path genome
    
    output:
    path("${genome.baseName}.mmi")
    
    when:
    params.isLong
    
    script:
    """
    minimap2 -ax splice -d ${genome.baseName}.mmi ${genome}
    """
}

//Mapping Nanopore reads to the genome and converting SAM file to coordinate-sorted BAM file
process MINIMAP2_MAP {
    container 'community.wave.seqera.io/library/minimap2_picard_samtools_star_subread:a1a7ac39cfcdefee'
    label 'multi_long'
    tag "on $read_id"
    publishDir "$params.outdir/Minimap2/", mode:'copy'
    
    input:
    tuple val(read_id), path(reads)
    path index
    
    output:
    path("map_${read_id}")
    tuple val(read_id), path("map_${read_id}/${read_id}_sorted.bam"), emit: sorted
    tuple val(read_id), path("map_${read_id}/${read_id}.bam.bai"), emit: index
    path("map_${read_id}/${read_id}.txt"), emit: stat
    
    when:
    params.isLong
        
    script:
    """
    mkdir -p map_${read_id} 
    paftools.js gff2bed ${params.annotation} > annotation.bed
    minimap2 -t {task.cpus} -ax splice ${index} --junc-bed annotation.bed ${reads} | samtools view -hbo ${read_id}.bam
    samtools sort -o map_${read_id}/${read_id}_sorted.bam ${read_id}.bam
    samtools index map_${read_id}/${read_id}_sorted.bam map_${read_id}/${read_id}.bam.bai
    samtools stats map_${read_id}/${read_id}_sorted.bam > map_${read_id}/${read_id}.txt
    """
}
