process DESEQ2 {
    conda "envs/r.yaml"
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