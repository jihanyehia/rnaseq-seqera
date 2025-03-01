//Run gene set enrichment analysis using fgsea
process FGSEA {
    debug true //to show echo in terminal
    container 'community.wave.seqera.io/library/bioconductor-biomart_bioconductor-deseq2_bioconductor-fgsea_r-ashr_pruned:b4d326d8b7640229'
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
