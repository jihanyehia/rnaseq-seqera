//Download fgsea gene sets of interest (species, collections, subcategories) 
process FGSEA_SETS {
    container 'community.wave.seqera.io/library/bioconductor-biomart_bioconductor-deseq2_bioconductor-fgsea_r-ashr_pruned:b4d326d8b7640229'
    publishDir "$params.outdir/FGSEA/", mode:'copy'
    
    output:
    path("gene_sets_dir"), emit: r_object
    path("gene_sets.csv"), emit: gene_sets
    
    script:
    """
    geneSets.R --collections "${params.collections}" --subcategory "${params.subcategories}" --species "${params.species}"
    """
}
