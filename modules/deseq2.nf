process DESEQ2 {
    container 'community.wave.seqera.io/library/bioconductor-biomart_bioconductor-deseq2_bioconductor-fgsea_r-ashr_pruned:b4d326d8b7640229'
    publishDir "$params.outdir/DESeq2/", mode:'copy'
    
    input:
    path(counts)
    
    output:
    path("normalized_counts.csv"), emit: norm_counts
    path("DESeq_results_*.csv"), emit: DE_genes
    path("DESeq_summary.txt")
    path("Downregulated_*.csv")
    path("Upregulated_*.csv")
    path("GSEA_metric_*.rnk"), emit: ranked
    path("GSEA_shrunk_*.rnk"), emit: ranked_shrunk
    path("GSEA_DE_*.rnk"), emit: ranked_DE
    path("Rplots.pdf")
    
    script:
    """
    DeSeq.R --counts ${counts} --metadata ${params.metadata} --pvalue ${params.pvalue} --lfcThreshold ${params.lfcThreshold}
    """
}
