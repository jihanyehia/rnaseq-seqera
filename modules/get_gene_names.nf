process GENE_NAMES {
    conda "envs/r.yaml"
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