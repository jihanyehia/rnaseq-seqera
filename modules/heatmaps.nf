//Create heatmaps of genes that were found to be differentially expressed in pathways that are among the top enriched
process HEATMAPS {
    debug true //to show echo in terminal
    conda "envs/py.yaml"
    label 'high'
    cpus '4'
    publishDir "$params.outdir/FGSEA/Heatmaps/", mode:'copy'
    
    input:
    path(counts)
    path(DE_genes)
    path(gene_sets)
    path(top_pathways)
    
    output:
    path("heatmaps.pdf"), optional: true
    
    script:
    """
    lines="\$(cat $DE_genes | wc -l)"
    if [[ \$lines -gt 1 ]]; then
        echo "DE genes detected. Heatmaps generated."
        python $projectDir/bin/heatmaps.py --counts ${counts} --g ${DE_genes} --gs ${gene_sets} --pathways ${top_pathways}
    else
        echo "No DE genes detected. Heatmaps omitted."
    fi
    """
}