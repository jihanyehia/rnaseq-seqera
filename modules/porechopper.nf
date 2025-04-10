//Adapter trimming for Nanopore
process PORECHOP {
    conda "envs/trim.yaml"
    label 'multi_long'
    tag "on $read_id"
    publishDir "$params.outdir/Porechop/", mode:'copy'
    
    input:
    tuple val(read_id), path(reads)
    
    output:
    tuple val(read_id), path("${read_id}_trimmed.fastq")
    
    when:
    params.isLong
    
    script:
    """
    porechop -t ${task.cpus} --discard_middle -i ${reads} -o ${read_id}_trimmed.fastq
    """
}

//Adapter filtering for Nanopore
process CHOPPER {
    conda "envs/trim.yaml"
    label 'multi_long'
    tag "on $read_id"
    publishDir "$params.outdir/CHOPPER/", mode:'copy'
    
    input:
    tuple val(read_id), path(trimmed)
    
    output:
    tuple val(read_id), path("${read_id}_chopped.fastq")
    
    when:
    params.isLong
    
    script:
    """
    cat ${trimmed} | chopper -q ${params.long_minQ} -l ${params.long_minLen} > ${read_id}_chopped.fastq
    """
}