//Generating MultiQC report
process MULTIQC {
    container 'community.wave.seqera.io/library/fastqc_multiqc_nanoplot:60ad77caa9e5f469'
    publishDir "$params.outdir/MultiQC/", mode:'copy'
    
    input:
    path multiqc_input_files
    val report_name
    
    output:
    path("${report_name}.html")
    
    script:
    """
    multiqc -n ${report_name}.html ${multiqc_input_files.join(' ')}
    """
}
