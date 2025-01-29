//Run gene set enrichment analysis using fgsea
process FGSEA {
    debug true //to show echo in terminal
    conda "envs/r.yaml"
    label 'high'
    cpus '4'
    publishDir "$params.outdir/FGSEA/", mode:'copy'
    
    input:
    path(ranks)
    path(gene_sets_dir)
    
    output:
    path("Rplots_*.pdf"), optional: true
    path("*_top_pathways_*.csv"), optional: true, emit: top_pathways
    
    script:
    """
    fgsea.R --file ${ranks} --gs ${gene_sets_dir} --nperm ${params.nperm}
    """
}