//Generating MultiQC report
process MULTIQC {
    conda "envs/qc.yaml"
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