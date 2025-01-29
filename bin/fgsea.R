#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))

pdf.options(width=18 , height=8)

#Function that runs fgsea on ranked lists of genes and generates a CSV file with the top 10 upregulated and downregulated pathways and PDF enrichment plots for each collection/subcategory of gene sets 
fgsea_results <- function(ranks, gene_set, nperm, coll_names, subcat_names, comparison, metric, source){
    #Group gene IDs by gene set names
    gene_list <- split(x = gene_set$gene_symbol, f = gene_set$gs_name)
    #Run fgsea
    res <- fgsea(gene_list, ranks, nPermSimple = nperm)
    #Get top 20 enriched pathways
    topUp <- res %>% filter(ES > 0) %>% arrange(padj) %>% head(n=10)
    topDown <- res %>% filter(ES < 0) %>% arrange(padj) %>% head(n=10)
    top <- bind_rows(topUp, topDown) %>% arrange(-NES)
    top_paths <- top[,-c(2:8)]
    #Generate file name and plot title
    if (unique(gene_set$gs_subcat) == "" | length(unique(gene_set$gs_subcat)) > 1) { # Collection without subcategories or no user-specified subcategory
        filename <- paste(as.character(gene_set$gs_cat[1]), "_top_pathways_", comparison, "_", source, ".csv", sep = "")
        title <- paste("Enrichment plot for the", as.character(coll_names[gene_set$gs_cat[1]]), "\n using ", metric, "\n", comparison)
    } else { # Collection with subcategory specified by user
        filename <- paste(as.character(gene_set$gs_cat[1]), "_", as.character(gene_set$gs_subcat[1]), "_top_pathways_", comparison, "_", source, ".csv", sep = "")
        title <- paste("Enrichment plot for the", as.character(coll_names[gene_set$gs_cat[1]]), as.character(subcat_names[gene_set$gs_subcat[1]]), "\n using ", metric, "\n", comparison)
    }
    write.csv(top_paths,file=filename,row.names=F)
    grid.arrange(left = title, plotGseaTable(gene_list[top$pathway], ranks, res, gseaParam=0.5, colwidths = c(5,3,0.8,0.8,0.8)))
}

option_list = list(optparse::make_option("--file",
                                        type = 'character',
                                        help = "path to ranked list file"))
option_list = c(option_list, list(optparse::make_option("--gs",
                                        type = 'character',
                                        help = "directory holding list of gene sets as R object")))
option_list = c(option_list, list(optparse::make_option("--nperm",
                                        type = 'integer',
                                        help = "Number of permutations for FGSEA")))

#Parse arguments
opt = optparse::parse_args(optparse::OptionParser(option_list = option_list))

ranks <- read.table(opt$file, header = TRUE, colClasses = c("character", "numeric"))
#assign gene IDs to metrics
ranks <- setNames(ranks$metric, ranks$Gene_ID)

gene_sets <- readRDS(file.path(opt$gs, 'gene_sets.rds'))

nperm <- opt$nperm


#Collection Names to be used in plot titles
coll_names <- c(
    "H"="Hallmark (H)",
    "C1"="Positional (C1)",
    "C2"="Curated (C2)",
    "C3"="Regulatory (C3)",
    "C4"="Computational (C4)",
    "C5"="Ontology (C5)",
    "C6"="Oncogenic (C6)",
    "C7"="Immunologic (C7)",
    "C8" = "Cell Type (C8)"
)

#Subcategory names to be used in plot titles
subcat_names <- c(
    #C2 subcategories
    "CGP" = "Chemical and Genetic Perturbations",
    "CP:BIOCARTA" = "Canonical Pathways from BioCarta Database",
    "CP:KEGG" = "Canonical Pathways from KEGG_MEDICUS Database",
    "CP" = "All Canonical Pathways",
    "CP:PID" = "Canonical Pathways from PID Database",
    "CP:REACTOME" = "Canonical Pathways from Reactome Database",
    "CP:WIKIPATHWAYS" = "Canonical Pathways from WikiPathways Database",
    #C3 subcategories
    "MIR:MIR_Legacy" = "miRNA Targets Legacy Set)",
    "TFT:TFT_Legacy" = "Transcription Factor Targets Legacy Set",
    "TFT:GTRD" = "Transcription Factor Targets from GTRD database",
    "MIR:MIRDB" = "miRNA targets from miRDB database",
    #C4 subcategories
    "CGN" = "Cancer Gene Neighborhoods", 
    "CM" = "Cancer Modules",
    #C5 subcategories
    "GO:BP" = "GO Biological Process", 
    "GO:CC" = "GO Cellular Component",
    "GO:MF" = "GO Molecular Function",
    "HPO" = "Human Phenotype Ontology",
    #C7 subcategories
    "VAX" = "Vaccine Response",
    "IMMUNESIGDB" = "Immunologic Signature from ImmuneSigDB"    
)

#Check which metric was used to rank genes (formula, shrunk LFC, or shrunk LFC on DE genes only)
filename = basename(opt$file)
if (grepl("metric", filename)) {
    comparison <- sub("GSEA_metric_(.*)\\.rnk", "\\1", filename)
    metric <- "-log10({p-adjusted}) * sign({Fold Change})"
    source <- "metric"
} else if (grepl("shrunk", filename)) {
    comparison <- sub("GSEA_shrunk_(.*)\\.rnk", "\\1", filename)
    metric <- "shrunk log2FC on all genes"
    source <- "shrunk"
} else if (grepl("DE", filename)) {
    comparison <- sub("GSEA_DE_(.*)\\.rnk", "\\1", filename)
    metric <- "shrunk log2FC on DE genes"
    source <- "DE"
} else {cat("Ranked file with unknown metric detected\n")}

pdf(file = paste0("Rplots_", comparison, "_", source, ".pdf"))

# Run analysis on each of the gene sets
if (is.data.frame(gene_sets)) {
    fgsea_results(ranks, gene_sets, nperm, coll_names, subcat_names, comparison, metric, source)
} else {
    for (gene_set in gene_sets) {
        fgsea_results(ranks, gene_set, nperm, coll_names, subcat_names, comparison, metric, source)
    }
}

#To save the results in a text format data:table::fwrite function can be used:
#fwrite(fgseaRes, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))