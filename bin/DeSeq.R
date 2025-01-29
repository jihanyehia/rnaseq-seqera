#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ashr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))

option_list = list(optparse::make_option("--counts",
                                        type = 'character',
                                        help = "path to count file"))

option_list = c(option_list, list(optparse::make_option("--metadata",
                                        type = 'character',
                                        help = "path to metadata file")))

option_list = c(option_list, list(optparse::make_option("--pvalue",
                                        type = 'double',
                                        help = "significance threshold for adjusted p-value")))

option_list = c(option_list, list(optparse::make_option("--lfcThreshold",
                                        type = 'double',
                                        help = "log2 fold change threshold")))

opt = optparse::parse_args(optparse::OptionParser(option_list = option_list))

#Read counts file with "Geneid" column as the row names
counts_table <- read.table(opt$counts, header = TRUE, sep = "\t", row.names = 1, quote = "") #to ignore quotes in file

#Remove first 5 columns (chr,start,end,strand,length) and _dedup.bam from sample names
counts_table <- counts_table[,-c(1:5)]
colnames(counts_table) <- sub("_dedup.bam", "", colnames(counts_table))
# Handle missing values by replacing NAs with zero
counts_table[is.na(counts_table)] <- 0

#Drop unwanted samples
#if(any(grepl("WE_A", colnames(counts_table)))){
#counts_table = subset(counts_table, select = -WE_A )}

#Read metadata file with "Sample name" column as the row names
metadata <- read.table(opt$metadata, header = TRUE, sep = "\t", row.names = 1)

#Extract reference levels, which is the last row of metadata dataframe, into a named list, then remove the row from the dataframe
reference_levels <- metadata[nrow(metadata), , drop = FALSE] #drop=FALSE ensures that the result remains a data.frame even if it has only one row or column
reference_list <- as.list(reference_levels)
metadata <- metadata[-nrow(metadata), , drop = FALSE]

# Reorder metadata rows to match the order of colnames in counts_table
metadata <- metadata[match(colnames(counts_table), rownames(metadata)), , drop = FALSE]

#Factor and relevel metadata to be used as colData
for (column in colnames(metadata)) {
    metadata[[column]] <- factor(metadata[[column]])
    ref_level <- reference_list[[column]]
    metadata[[column]] <- relevel(metadata[[column]], ref = ref_level)
}

#Build dataset, create contrasts,and run differential expression  
if (ncol(metadata) == 1) {
    design_formula <- as.formula(paste("~", colnames(metadata)))
    dds <- DESeqDataSetFromMatrix(countData = counts_table, colData = metadata, design = design_formula)
    cond_levels <- levels(metadata[[colnames(metadata)]])
    #Pairwise comparisons with reference group (no comparison between non-reference groups)
    contrasts <- list()
    for (i in 2:length(cond_levels)) {
        # Create contrast for each level against the reference (cond_levels[1])
        contrasts[[i - 1]] <- c(colnames(metadata), cond_levels[i], cond_levels[1])
    }
    contrast_names <- paste(cond_levels[2:length(cond_levels)], "vs", cond_levels[1], sep = "_")
} else if (ncol(metadata) > 1) {
    #If there is more than one factor under study, combine factors together to build dataset based on these combinations
    #Supports only two conditions per factor
    combined_factors <- apply(metadata, 1, function(row) paste(row, collapse = "_"))
    metadata$combined <- factor(combined_factors) 
    dds <- DESeqDataSetFromMatrix(countData = counts_table, colData = metadata, design = ~ combined)
    columns <- colnames(metadata)
    cond1_levels <- levels(metadata[[columns[1]]])
    cond2_levels <- levels(metadata[[columns[2]]])
    contrasts <- list(
        c("combined", paste(cond1_levels[2], cond2_levels[1], sep = "_"), paste(cond1_levels[1], cond2_levels[1], sep = "_")),       # Compring cond1 levels with cond2's reference level
        c("combined", paste(cond1_levels[2], cond2_levels[2], sep = "_"), paste(cond1_levels[1], cond2_levels[2], sep = "_")),   # Compring cond1 levels with cond2's non-reference level
        c("combined", paste(cond1_levels[1], cond2_levels[2], sep = "_"), paste(cond1_levels[1], cond2_levels[1], sep = "_")),     # Compring cond2 levels with cond1's reference level
        c("combined", paste(cond1_levels[2], cond2_levels[2], sep = "_"), paste(cond1_levels[2], cond2_levels[1], sep = "_"))     # Compring cond2 levels with cond1's non-reference level
    )
    contrast_names <- c(
    paste(cond1_levels[2], "vs", cond1_levels[1], cond2_levels[1], sep = "_"),  # e.g., BPA_vs_ctl_WT
    paste(cond1_levels[2], "vs", cond1_levels[1], cond2_levels[2], sep = "_"),  # e.g., BPA_vs_ctl_Fmr1
    paste(cond2_levels[2], "vs", cond2_levels[1], cond1_levels[1], sep = "_"),  # e.g., Fmr1_vs_WT_ctl
    paste(cond2_levels[2], "vs", cond2_levels[1], cond1_levels[2], sep = "_")   # e.g., Fmr1_vs_WT_BPA
    )
}
names(contrasts) <- contrast_names
dds <- DESeq(dds)

pval <- opt$pvalue
lfcThreshold <- opt$lfcThreshold
                            
#Loop over contrasts, extract results with an adjusted p-value cut-off set by user (default=0.05), and store results in corresponding files
for (i in seq_along(contrasts)){
    cont <- contrasts[[i]]
    res <- results(dds, contrast = cont, alpha = pval)  
    res_shrunk <- lfcShrink(dds, res = res, type="ashr")
    
    #Filter DE genes
    DE_genes <- as.data.frame(res) %>%
    filter(abs(log2FoldChange) >= lfcThreshold, padj<=pval) %>%
    arrange(padj)
    
    #Upregulated genes
    up_genes <- DE_genes %>%
    filter(log2FoldChange > lfcThreshold) %>%
    arrange(padj, desc(log2FoldChange))
    
    #Downregulated genes
    down_genes <- DE_genes %>%
    filter(log2FoldChange < lfcThreshold) %>%
    arrange(padj, log2FoldChange)
    
    #Prepare ranked list of genes sorted by the metric formula "-log10({p value}) * sign({Fold Change})" for gsea
    df <- as.data.frame(res)
    df <- df[!is.na(df$pvalue),] #get rid of NA in p-value
    df <- df[!is.na(df$padj),] #get rid of NA in padj
    df$pvalue[df$padj == 0] <- 1e-50  #change p-adj of 0 to a very low number to prevent Inf in metrics
    df$signFC <- sign(df$log2FoldChange)
    df$logP <- -log10(df$padj)
    df$metric <- df$logP * df$signFC
    ranked_metric <- data.frame(row.names(df), df$metric)
    colnames(ranked_metric) = c("Gene_ID", "metric")
    ranked_metric <- ranked_metric[order(-ranked_metric$metric), ]
    
    #Prepare ranked list of shrunken LFC values for all genes for gsea
    df_shrunk <- as.data.frame(res_shrunk)
    df_shrunk <- df_shrunk[!is.na(df_shrunk$pvalue),] #get rid of NA in p-value
    df_shrunk <- df_shrunk[!is.na(df_shrunk$padj),] #get rid of NA in padj
    shrunk_gsea <- data.frame(rownames(df_shrunk), df_shrunk$log2FoldChange)
    colnames(shrunk_gsea) = c("Gene_ID", "metric")
    shrunk_gsea <- shrunk_gsea[order(-shrunk_gsea$metric), ]
    
    #Prepare ranked list of shrunken LFC values for DE genes only
    DE_shrunk <- df_shrunk %>%
    filter(abs(log2FoldChange) >= lfcThreshold, padj<=pval) %>%
    arrange(desc(log2FoldChange))
    shrunk_DE <- data.frame(rownames(DE_shrunk), DE_shrunk$log2FoldChange)
    colnames(shrunk_DE) = c("Gene_ID", "metric")
    
    write.table(DE_genes, file=paste0("DESeq_results_", contrast_names[i], ".csv"), sep=",", quote=F, col.names=NA)
    write.table(up_genes, file=paste0("Upregulated_", contrast_names[i], ".csv"), sep=",", quote=F, col.names=NA)
    write.table(down_genes, file=paste0("Downregulated_", contrast_names[i], ".csv"), sep=",", quote=F, col.names=NA)
    write.table(ranked_metric, file=paste0("GSEA_metric_", contrast_names[i], ".rnk"), sep="\t", quote=F, row.names=FALSE)
    write.table(shrunk_gsea, file=paste0("GSEA_shrunk_", contrast_names[i], ".rnk"), sep="\t", quote=F, row.names=FALSE)
    write.table(shrunk_DE, file=paste0("GSEA_DE_", contrast_names[i], ".rnk"), sep="\t", quote=F, row.names=FALSE)
  
    #Get condition to use it in plot titles and summary file
    info <- mcols(res)
    log2fc_description <- info[rownames(info) == "log2FoldChange", "description"]
    if (ncol(metadata) > 1){
        cond <- gsub(".* \\(MLE\\): combined ", "", log2fc_description)
    } else {
        cond <- gsub(paste0(".* \\(MLE\\): ", colnames(metadata), " "), "", log2fc_description) 
    }
    
    #Write results to output file
    sink(file="DESeq_summary.txt", append = TRUE)
    print(cond)
    summary(res)
    sink()
    
    #MA plot
    plotMA(res, ylim=c(-8,8), main=paste0("MA plot with no shrinkage for ", cond))
    plotMA(res_shrunk, ylim=c(-8,8), main=paste0("MA plot with shrinkage for ", cond))
  
    #P-value distribution
    hist(res$pvalue,main=paste0("P-value distribution for ", cond),xlab="P-values")
  
    #P-adjusted distribution
    hist(res$padj,main=paste0("P-adjusted distribution for ", cond),xlab="P-adjusted")
}
                
#Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts <- cbind(rownames(normalized_counts), normalized_counts) #to make the index the first column of dataframe
colnames(normalized_counts)[1] <- "Gene" #to rename the first column
                
write.table(normalized_counts, file="normalized_counts.csv", sep=",", quote=F, row.names=FALSE)          

                
#Transform dataset
rld <- rlog(dds, blind = TRUE) #blind=true: transformation unbiased to sample condition information

###Visualization###
#PCA plot
if (ncol(metadata) > 1){
    p <- plotPCA(rld, intgroup=c(colnames(metadata)[1], colnames(metadata)[2])) #intgroup should specify columns of colData(dds)
    print(p + geom_text_repel(aes_string(x = "PC1", y = "PC2", label = factor(colnames(dds))), size = 3)  + geom_point() + guides(color=guide_legend(title=colnames(metadata)[1]), shape=guide_legend(title=colnames(metadata)[2])))
} else{
    p <- plotPCA(rld, intgroup=colnames(metadata)) #intgroup should specify columns of colData(dds)
    print(p + geom_text_repel(aes_string(x = "PC1", y = "PC2", label = factor(colnames(dds))), size = 3)  + geom_point()) #requires ggplot2 and ggrepel
}


#Sample-sample comparison heatmap
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color=heat.colors, main = "Sample-to-Sample correlation heatmap",fontsize = 10, fontsize_row = 10, height = 20)

#Dispersion estimates
plotDispEsts(dds, main = "Dispersion estimates") #on the DESeq dataset