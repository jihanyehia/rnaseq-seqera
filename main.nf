#!/usr/bin/env nextflow

//use Nextflow DSL 2
nextflow.enable.dsl = 2
      
def helpMessage() {
  log.info """
   Usage:
    The typical command for running the pipeline is as follows:
        nextflow run main.nf

   Optional arguments:
         --project                Project name to be used as the name of the directory holding the output files
         --metadata               Full path to metadata file with samples in rows and conditions in columns. Last row (Reference) should contain the reference levels of each condition (e.g. Control, WT)
         --isLong                 Boolean to run pipeline in long read mode. Default is short read mode [false]
         --isSingle               Boolean to run pipeline in single-end read mode. Default is paired-end read mode [false]
         --PEreads                Full path to paired-end read files with Rex (example:*_{1,2}.fq.gz)
         --SEreads                Full path to single-end read files with Rex (example:*.fastq)
         --genome                 Reference genome (full path required)
         --annotation             Gene annotation file (full path required)
         --refFlat                Gene annotations in refFlat form (full path required)
         --outdir                 Parent directory to place intermediate and final output files
         --[short/long]_minQ      Minimum quality score for filtering short/long reads
         --[short/long]_minLen    Minimum read length for short/long filtering
         --overhang               STAR's sjdboverhang (should be max(read length) - 1)
         --strandedness           FeatureCounts check for performing strand-specific read counts (unstranded = 0, stranded = 1, reversely stranded = 2). Default = 0. 
         --species                Speices under study for MSigDB (Homo sapiens, Drosophila melanogaster, Mus musculus)
         --collections            MSigDB collections for pathway enrichment analysis (C1,C2,C3,C4,C5,C6,C7,C8,H)
         --subcategories          MSigDB subcategories of collections (used if interested in certain subcategories not all)
         --pvalue                 Significance threshold for adjusted p-values to identify differentially expressed genes. Default = 0.05
         --lfcThreshold           Log2 fold change threshold to filter genes based on expression changes. Default = 0
         --nperm                  FGSEA's number of permutations. Default = 1000
         --help                   This usage statement
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

log.info """\
             R N A S E Q  P I P E L I N E    
         ====================================
         Results: ${params.outdir}
         Reads : ${params.isSingle ? params.SEreads : params.PEreads}
         Genome: ${params.genome}
         Annotation: ${params.annotation}
         Long read mode: ${params.isLong}
         Single read mode: ${params.isSingle}
         """
         .stripIndent()
 
//Import modules
include { NANOPLOT ; NANOPLOT as NANOPLOT_FILTERED } from './modules/nanoplot.nf'
include { FASTQC ; FASTQC as FASTQC_FILTERED } from './modules/fastqc.nf'
include { MULTIQC ; MULTIQC as MULTIQC_FILTERED ; MULTIQC as MULTIQC_ALIGNED ; MULTIQC as MULTIQC_DUPS ; MULTIQC as MULTIQC_COUNTS } from './modules/multiqc.nf'
include { FASTP_SINGLE ; FASTP_PAIRED } from './modules/fastp.nf'
include { PORECHOP ; CHOPPER } from './modules/porechopper.nf'
include { MINIMAP2_INDEX ; MINIMAP2_MAP } from './modules/minimap2.nf'
include { STAR_INDEX ; STAR_MAP } from './modules/star.nf'
include { DUPLICATES } from './modules/duplicates.nf'
include { FEATURECOUNTS } from './modules/counts.nf'
include { GENE_NAMES } from './modules/get_gene_names.nf'
include { DESEQ2 } from './modules/deseq2.nf'
include { FGSEA_SETS } from './modules/fgsea_sets.nf'
include { FGSEA } from './modules/fgsea.nf'
include { HEATMAPS } from './modules/heatmaps.nf'


//Main workflow
workflow {
    def reads_ch
    if (params.isSingle) {
        //Convert the input from path into tuple(read_id, path)
        reads_ch = Channel
            .fromPath(params.SEreads, checkIfExists: true)
            .map { tuple( it.getBaseName(it.name.endsWith('.gz')? 2: 1), it ) }
    } else {
        reads_ch = Channel.fromFilePairs(params.PEreads, checkIfExists: true)
    }
    NANOPLOT(reads_ch, "raw")
    FASTQC(reads_ch, "raw")
    MULTIQC(FASTQC.out.collect(), "MultiQC_raw")
    def trimmed = (params.isSingle ? FASTP_SINGLE(reads_ch) : FASTP_PAIRED(reads_ch))
    PORECHOP(reads_ch)
    CHOPPER(PORECHOP.out)
    NANOPLOT_FILTERED(CHOPPER.out, "filtered")
    FASTQC_FILTERED(CHOPPER.out, "filtered")
    MULTIQC_FILTERED((params.isLong ? FASTQC_FILTERED.out.collect() : trimmed.json.collect()), "MultiQC_filtered")
    MINIMAP2_INDEX(params.genome)
    MINIMAP2_MAP(CHOPPER.out, MINIMAP2_INDEX.out)
    STAR_INDEX()
    STAR_MAP(STAR_INDEX.out, trimmed.filtered)
    MULTIQC_ALIGNED((params.isLong ? MINIMAP2_MAP.out.stat.collect() : STAR_MAP.out.log_final.collect()), "MultiQC_aligned")
    DUPLICATES(params.isLong ? MINIMAP2_MAP.out.sorted : STAR_MAP.out.sorted_bam)
    MULTIQC_DUPS(DUPLICATES.out.premetrics.mix(DUPLICATES.out.postmetrics).mix(DUPLICATES.out.dedupmetrics).collect(), "MultiQC_dups")
    FEATURECOUNTS(DUPLICATES.out.dedup_bam.collect())
    MULTIQC_COUNTS(FEATURECOUNTS.out.logs, "MultiQC_counts")
    //Pass the counts table to GENE_NAMES for gene ID conversion if the species under study is drosophila or mouse
    GENE_NAMES(FEATURECOUNTS.out.counts)
    DESEQ2(GENE_NAMES.out)
    FGSEA_SETS()
    combined_ranks = DESEQ2.out.ranked.mix(DESEQ2.out.ranked_shrunk).mix(DESEQ2.out.ranked_DE)
    FGSEA(combined_ranks.flatten(), FGSEA_SETS.out.r_object)
    //HEATMAPS(DESEQ2.out.norm_counts, DESEQ2.out.DE_genes.flatten(), FGSEA_SETS.out.gene_sets, FGSEA.out.top_pathways.collect())
}


workflow.onComplete { 
    log.info ( workflow.success ? "\nSuccessfully Completed!" : "Oops .. Something went wrong: ${workflow.errorMessage}" )
}
