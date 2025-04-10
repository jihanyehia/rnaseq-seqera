//Download fgsea gene sets of interest (species, collections, subcategories) 
process FGSEA_SETS {
    conda "envs/r.yaml"
    publishDir "$params.outdir/FGSEA/", mode:'copy'
    
    output:
    path("gene_sets_dir"), emit: r_object
    path("gene_sets.csv"), emit: gene_sets
    
    script:
    """
    geneSets.R --collections "${params.collections}" --subcategory "${params.subcategories}" --species "${params.species}"
    """
}