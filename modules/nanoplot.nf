//Quality control on Nanopore data (stage = raw | filtered)
process NANOPLOT {
    container 'community.wave.seqera.io/library/fastqc_multiqc_nanoplot:60ad77caa9e5f469'
    tag "on $read_id"
    publishDir "$params.outdir/Nanoplot/", mode:'copy'
    
    input:
    tuple val(read_id), path(reads)
    val stage
    
    output:
    path("${read_id}_${stage}")
    
    when:
    params.isLong
    
    script:
    """
    mkdir -p ${read_id}_${stage}
    NanoPlot --fastq ${reads} -o ${read_id}_${stage}
    """
}
