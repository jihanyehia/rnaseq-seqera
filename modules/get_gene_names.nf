process GENE_NAMES {
    container 'community.wave.seqera.io/library/bioconductor-biomart_bioconductor-deseq2_bioconductor-fgsea_r-ashr_pruned:b4d326d8b7640229'
    publishDir "$params.outdir/DESeq2/", mode:'copy'
    
    input:
    path(counts)
    
    output:
    path("updated_counts.txt")
    
    
    script:
    """
    if [[ "$params.species" == 'Drosophila melanogaster' || "$params.species" == 'Mus musculus' ]]; then
        getGeneNames.R --file ${counts} --species "${params.species}"
    else
        # Copy or rename the original counts if no gene name update is needed
        cp ${counts} updated_counts.txt
    fi
    """
}
