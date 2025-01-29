#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(biomaRt))

map_gene_names <- function (counts_table, info) {
    # Change index from gene IDs to gene names with biomaRt
    
    # Extract gene IDs from counts table and Ensembl
    count_gene_ids <- rownames(counts_table)
    ensmbl_gene_id <- info$attributes[1]
    # Initialize the Ensembl mart and dataset
    ensembl <- useMart("ensembl")
    species <- useDataset(info$dataset, mart = ensembl)
    # Query to map gene IDs to gene names
    result <- getBM(
        attributes = info$attributes,
        filters = ensmbl_gene_id,
        values = count_gene_ids,
        mart = species
    )
    
    # Remove duplicates from result
    result <- result[!duplicated(result[[ensmbl_gene_id]]), ]

    # Make sure order matches that in the counts table
    result_ordered <- result[match(count_gene_ids, result[[ensmbl_gene_id]]), ]
    
    # Replace NAs in external_gene_name with unique identifiers (NA_1, NA_2...)
    result_ordered$external_gene_name[is.na(result_ordered$external_gene_name)] <- paste0("NA_", seq_along(which(is.na(result_ordered$external_gene_name))))
    
    # Add gene_name to the counts table to be set as index and ensure uniqueness
    counts_table$gene_name <- make.unique(result_ordered$external_gene_name)

    # Set gene_name column as index and remove column
    rownames(counts_table) <- counts_table$gene_name
    counts_table <- counts_table [,!(names(counts_table) %in% "gene_name")]
    
    return(counts_table)
}

# Define a lookup table for species and their information: species = (dataset, attributes = c(filter, gene_name))
species_lookup <- list(
  "Drosophila melanogaster" = list(dataset = "dmelanogaster_gene_ensembl", attributes = c("flybase_gene_id", "external_gene_name")),
  "Mus musculus" = list(dataset = "mmusculus_gene_ensembl", attributes = c("ensembl_gene_id", "external_gene_name"))
)

option_list = list(optparse::make_option("--file",
                                        type = 'character',
                                        help = "path to count file"))

option_list = c(option_list, list(optparse::make_option("--species",
                                        type = 'character',
                                        help = "Species under study")))

opt = optparse::parse_args(optparse::OptionParser(option_list = option_list))

#Read counts file with "Geneid" column as the row names
counts_table <- read.table(opt$file, header = TRUE, sep = "\t", row.names = 1)

species <- opt$species
# Handle cases where the species is not in the lookup table
if (!(species %in% names(species_lookup))) {
    stop("Species not supported or not found in lookup table: ", species)
}
# Retrieve the info based on the species
info <- species_lookup[[species]]

# Handle cases where the species is not in the lookup table
if (!(species %in% names(species_lookup))) {
    stop("Species not supported or not found in lookup table: ", species)
}

updated_counts_table <- map_gene_names(counts_table, info)

#Save the updated counts table as TXT file to be read by DeSeq2.R 
write.table(updated_counts_table, file = "updated_counts.txt", sep = "\t", row.names = TRUE, quote = FALSE)