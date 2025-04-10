//Adapter trimming and filtering for short single-end reads (default min_quality_score=15; default min_length=15; low complexity filter; threading;json report)
process FASTP_SINGLE {
    conda "envs/trim.yaml"
    label 'multi_short'
    clusterOptions "--nodes=1"
    tag "on $read_id"
    publishDir "$params.outdir/Fastp/", mode:'copy', pattern: '*trimmed.fq.gz'
    
    input:
    tuple val(read_id), path(reads)
    
    output:
    tuple val(read_id), path("${read_id}_trimmed.fq.gz"), emit: filtered
    path("${read_id}.fastp.json"), emit: json
    
    when:
    !params.isLong
    
    script:
    """
    fastp -i ${reads[0]} -o ${read_id}_trimmed.fq.gz -q ${params.short_minQ} -l ${params.short_minLen} -y --thread ${task.cpus} --json ${read_id}.fastp.json
    """

}

//Adapter trimming and filtering for short paired-end reads (default min_quality_score=15; default min_length=15; low complexity filter; threading; automatically detect adapters for paired-end reads; json report)
process FASTP_PAIRED {
    conda "envs/trim.yaml"
    label 'multi_short'
    clusterOptions "--nodes=1"
    tag "on $read_id"
    publishDir "$params.outdir/Fastp/", mode:'copy', pattern: '*trimmed.fq.gz'
    
    input:
    tuple val(read_id), path(reads)
    
    output:
    tuple val(read_id), path("${read_id}_{1,2}_trimmed.fq.gz"), emit: filtered
    path("${read_id}.fastp.json"), emit: json
    
    when:
    !params.isLong
    
    script:
    """
    fastp -i ${reads[0]} -o ${read_id}_1_trimmed.fq.gz -I ${reads[1]} -O ${read_id}_2_trimmed.fq.gz -q ${params.short_minQ} -l ${params.short_minLen} -y --thread ${task.cpus} --detect_adapter_for_pe --json ${read_id}.fastp.json
    """

}