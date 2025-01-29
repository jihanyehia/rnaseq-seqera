//Dealing with duplicates: AlignmentSummaryMetrics, MarkDuplicates, AlignmentSummaryMetrics
process DUPLICATES {
    conda "envs/counts.yaml"
    tag "on $read_id"
    label params.isLong ? 'multi_long' : 'high'  
    publishDir "$params.outdir/Duplicates/", mode:'copy'
    
    input:
    tuple val(read_id), path(bam)
   
   output:
    path("${read_id}_dup/${read_id}_metrics_before"), emit: premetrics
    path("${read_id}_dup/${read_id}_dedup.bam"), emit: dedup_bam
    path("${read_id}_dup/${read_id}.marked_dup_metrics.txt"), emit: dedupmetrics
    path("${read_id}_dup/${read_id}_metrics_after"), emit: postmetrics
    
    script:
    """
    mkdir -p ${read_id}_dup/${read_id}_metrics_before ${read_id}_dup/${read_id}_metrics_after
    #Check if read group line exists in bam files and add placeholders if it doesn't
    if picard ViewSam I=${bam} | grep -q '^@RG'; then
        cp ${bam} ${read_id}_dup/${read_id}_with_rg.bam
    else
        picard AddOrReplaceReadGroups \
            I=${bam} \
            O=${read_id}_dup/${read_id}_with_rg.bam \
            RGID=${read_id} \
            RGLB=lib1 \
            RGPL=UNKNOWN \
            RGPU=unit1 \
            RGSM=${read_id} \
            CREATE_INDEX=true
    fi
    #Collect pre-metrics
    picard -Xmx10G CollectMultipleMetrics \
        I=${read_id}_dup/${read_id}_with_rg.bam \
        O=${read_id}_dup/${read_id}_metrics_before/${read_id}_metrics_before \
        R=${params.genome} \
        PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=RnaSeqMetrics \
        REF_FLAT=${params.refFlat}
    #Mark and remove duplicates
    picard -Xmx10G MarkDuplicates \
        I=${read_id}_dup/${read_id}_with_rg.bam \
        O=${read_id}_dup/${read_id}_dedup.bam \
        M=${read_id}_dup/${read_id}.marked_dup_metrics.txt \
        REMOVE_DUPLICATES=true \
        ASSUME_SORT_ORDER=coordinate \
        CREATE_INDEX=true
    #Collect post-metrics
    picard -Xmx10G CollectMultipleMetrics \
        I=${read_id}_dup/${read_id}_dedup.bam \
        O=${read_id}_dup/${read_id}_metrics_after/${read_id}_metrics_after \
        R=${params.genome} \
        PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=RnaSeqMetrics \
        REF_FLAT=${params.refFlat}
    """
}