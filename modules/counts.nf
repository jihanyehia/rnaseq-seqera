//Quantifying reads ( -L long-read mode, -C chimeric reads are not counted, -p reads are paired, -B only reads with both ends successfully mapped are counted)
process FEATURECOUNTS {
    container 'community.wave.seqera.io/library/minimap2_picard_samtools_star_subread:a1a7ac39cfcdefee'
    label params.isLong ? 'multi_long' : 'multi_short'
    publishDir "$params.outdir/FeatureCounts/", mode:'copy'
    
    input:
    path(bams)
    
    output:
    path("counts.txt"), emit: counts
    path("counts.txt.summary"), emit: logs
    
    script:
    """
    featureCounts ${params.isSingle ? (params.isLong ? '-L' : '') : '-p'} \
    -C -T 32 -t exon -g gene_id -a ${params.annotation} -s ${params.strandedness} -o counts.txt ${bams}
    """
}
