//Quality control (stage = raw | filtered)
process FASTQC {
    conda "${projectDir}/envs/qc.yaml"
    tag "on $read_id"
    label params.isLong ? 'multi_long' : 'multi_short'
    publishDir "$params.outdir/FastQC/", mode:'copy'
    
    input:
    tuple val(read_id), path(reads)
    val stage
    
    output:
    path("${read_id}_${stage}")
    
    script:
    """
    mkdir -p ${read_id}_${stage}
    fastqc -t ${task.cpus} -o ${read_id}_${stage} -q ${reads}
    """
}
