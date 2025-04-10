#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(biomaRt)) # map gene IDs to gene names using Ensembl
suppressPackageStartupMessages(library(dplyr))


aggregate_counts_by_gene <- function(counts_df, gene_map,
                                     id_col = "ensembl_gene_id",
                                     name_col = "external_gene_name") {

  
  # Identify non-sample metadata columns (excluding ensembl_id) from sample columns
  sample_cols <- grep("_dedup.bam", colnames(counts_df), value = TRUE)
  meta_cols <- setdiff(colnames(counts_df), c(sample_cols, "ensembl_id"))

  # Join counts data with gene mapping (biomaRt results) by gene ID.
  # Filter out rows where the gene name is missing (NA).
  joined_df <- counts_df %>%
    left_join(gene_map, by = c("ensembl_id" = id_col)) %>%
    filter(!is.na(!!sym(name_col)) & !!sym(name_col) != "")

  # For each gene name, record which Ensembl IDs contributed to it.
  # This helps with traceability in case multiple Ensembl IDs map to the same gene name.
  id_trace <- joined_df %>%
    group_by(!!sym(name_col)) %>%
    summarise(ensembl_ids = paste(ensembl_id, collapse = ";"), .groups = "drop")

  # Keep metadata from the first entry per gene name
  meta_first <- joined_df %>%
    group_by(!!sym(name_col)) %>%
    summarise(across(all_of(meta_cols), ~ first(.x)), .groups = "drop")

  # Aggregate numeric count columns by gene name.
  # This collapses rows with the same gene name by summing counts across Ensembl IDs.
  counts_agg <- joined_df %>%
    group_by(!!sym(name_col)) %>%
    summarise(across(all_of(sample_cols), sum), .groups = "drop")

  # Merge all into final result: metadata, sample counts, ensembl_id
  result <- meta_first %>%
    left_join(counts_agg, by = name_col) %>%
    left_join(id_trace, by = name_col)

  # Set rownames as gene names, after converting tibble to dataframe (allows assigning rownames)
  result <- as.data.frame(result)
  rownames(result) <- result[[name_col]]
  result[[name_col]] <- NULL

  return(result)
}

    
map_gene_names <- function (counts_table, info) {
    # Change index from gene IDs to gene names with biomaRt
    
    # Extract gene IDs from counts table
    count_gene_ids <- rownames(counts_table)
    ensmbl_gene_id <- info$attributes[1]
    # Initialize biomaRt
    ensembl <- useMart("ensembl")
    species <- useDataset(info$dataset, mart = ensembl)
    # Query to map gene IDs to gene names
    gene_map <- getBM(
        attributes = info$attributes,
        filters = ensmbl_gene_id,
        values = count_gene_ids,
        mart = species
    )
    
    # Aggregate using the wrapper
    counts_table$ensembl_id <- rownames(counts_table)
    aggregated <- aggregate_counts_by_gene(counts_table, gene_map,
                                         id_col = ensmbl_gene_id,
                                         name_col = info$attributes[2])
    
    return(aggregated)
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

# Read counts file with "Geneid" column as the row names
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

# Save the updated counts table as TXT file to be read by DeSeq2.R 
write.table(updated_counts_table, file = "updated_counts.txt", sep = "\t", row.names = TRUE, quote = FALSE)