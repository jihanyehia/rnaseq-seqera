//Inferring strandedness and quantifying reads ( -L long-read mode, -C chimeric reads are not counted, -p reads are paired, -B only reads with both ends successfully mapped are counted)
process FEATURECOUNTS {
    debug true //to show echo in terminal
    conda "envs/counts.yaml"
    label params.isLong ? 'multi_long' : 'multi_short'
    publishDir "$params.outdir/FeatureCounts/", mode:'copy'
    
    input:
    path(bams)
    
    output:
    path("counts.txt"), emit: counts
    path("counts.txt.summary"), emit: logs
    
    script:
    """
    #Check if string length is zero
    if [[ -z "${params.strandedness}" ]]; then
        echo "Strandedness not provided, running infer_experiment.py..."
        
        #Convert annotation to bed format for strandedness
        gtf2bed < ${params.annotation} > annotation.bed

        #Use the first BAM file for inferring strandedness
        infer_experiment.py -r annotation.bed -i $bams[0] > strandedness.txt

        #Parse inferred strandedness
        strand=\$(awk '
        NR==4 {s1=\$NF}  # Extract the fraction from the second line (after colon)
        NR==5 {s2=\$NF}  # Extract the fraction from the third line (after colon)
        END {
            if (s1 > 0.9) print 1;  # Forward stranded
            else if (s2 > 0.9) print 2;  # Reverse stranded
            else print 0;  # Unstranded
        }' strandedness.txt)

        echo "Detected strandedness: \$strand"

    else
        strand="${params.strandedness}"  # Use the provided strandedness
    fi

    #Run featurecounts
    featureCounts ${params.isSingle ? (params.isLong ? '-L' : '') : '-p'} \
    -C -T 32 -t exon -g gene_id -a ${params.annotation} -s \$strand -o counts.txt ${bams}
    """
}