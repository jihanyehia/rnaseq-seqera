#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(msigdbr)) #retrieve gene sets from MSigDB database
suppressPackageStartupMessages(library(dplyr)) #manipulate data(grouping and summarizing gene set data)

# Function to find the category for a given subcategory code
find_category <- function(subcat, subcategory_dict) {
    for (collection in names(subcategory_dict)) {
        if (subcat %in% subcategory_dict[[collection]]) {
            return(collection)
        }
    }
    return(NULL)  # Return NULL if category not found
}

#Create arguments that the user can use to specify their gene sets and species of interest
option_list <- list(optparse::make_option("--collections", 
  type = "character", 
  help = "List of collections separated by commas"
))
option_list = c(option_list, list(optparse::make_option("--subcategory",
                                        type = 'character',
                                        help = "List of subcategories separated by commas")))
option_list = c(option_list, list(optparse::make_option("--species",
                                        type = 'character',
                                        help = "Name of species under study")))

# Parse the argument
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
if (!is.null(args$collections)) {
    if (grepl(",", args$collections)) {
        collections <- unlist(strsplit(args$collections, ","))
        collections <- trimws(collections)
    } else {
        collections <- args$collections
    }
    collections <- as.vector(collections)
} else {
    stop("At least one collection must be provided.")
}
if (!is.null(args$subcategory)) {
    if (grepl(",", args$subcategory)) {
        subcats <- unlist(strsplit(args$subcategory, ","))
        subcats <- trimws(subcats)
    } else {
        subcats <- args$subcategory
        subcats <- trimws(subcats)
    }
    subcats <- as.vector(subcats)
} else {
    subcats <- list() #Create an empty list in case no subcategory is provided
}
species <- args$species

# Create a dictionary-like structure with subcategories grouped by key
subcategory_dict <- list(
    C2 = c("CGP", "CP:BIOCARTA", "CP:KEGG", "CP", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS"),
    C3 = c("MIR:MIR_Legacy", "TFT:TFT_Legacy", "TFT:GTRD", "MIR:MIRDB"),
    C4 = c("CGN", "CM"),
    C5 = c("GO:BP", "GO:CC", "GO:MF", "HPO"),
    C7 = c("VAX", "IMMUNESIGDB")
)

# Get gene sets and combine them in one list
if (length(collections) == 1 ) {
    if (length(subcats) == 0) {
        gene_sets <- msigdbr(species = species, category = collections[1])
    } else if (length(subcats) == 1) {
        gene_sets <- msigdbr(species = species, category = collections[1], subcategory = subcats[1])
    } else { # Multiple subcategories 
        gene_sets <- lapply(subcats, function(sub) msigdbr(species = species, category = collections[1], subcategory = sub))
    }
} else { # Multiple Collections
    if (length(subcats) == 0) { # No subcategories specified
        gene_sets <- lapply(collections, function(col) msigdbr(species = species, category = col))
    } else { # One or more subcategories specified
       coll_handled <- c() # To keep track of collections with specified subcategories
        gene_sets <- lapply(subcats, function(subcat) {
            coll <- find_category(subcat, subcategory_dict)
            coll_handled <<- unique(c(coll_handled, coll)) # Use <<- to modify the outer scope variable
            msigdbr(species = species, category = coll, subcategory = subcat)})
        unhandled_coll <- setdiff(collections, coll_handled)
        if (length(unhandled_coll) > 0) { # Run msigdbr on remaining collections without specified subcategories
            gene_sets_unhandled <- lapply(unhandled_coll, function(col) {
                msigdbr(species = species, category = col)})
            gene_sets <- c(gene_sets, gene_sets_unhandled)
        } 
    }
}
    
#Define the directory where the objects to be used in fgsea.R will be saved and save the objects to RDS files in the specified directory
directory <- "gene_sets_dir"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
saveRDS(gene_sets, file.path(directory, 'gene_sets.rds'))
 

#Build a table dataframe of all collections with category (gs_cat), subcategory (gs_subcat), pathway (gs_name), and a list of its respective genes (gene_symbols) to be used in the python script for heatmaps 
table_all <- data.frame(gs_cat = character(), gs_subcat = character(), gs_name = character(), gene_list = I(list()))
if (is.data.frame(gene_sets)) { #only one gene set
  table_all <- gene_sets %>%
    group_by(gs_name, gs_cat, gs_subcat) %>%
    summarise(gene_list = list(unique(gene_symbol))) %>%
    ungroup() %>% select(gs_cat, gs_subcat, gs_name, gene_list)
  table_all$gene_list <- sapply(table_all$gene_list, toString)
} else { #more than one gene set
  for (gene_set in gene_sets) {
    table <- gene_set %>%
      group_by(gs_name, gs_cat, gs_subcat) %>%
      summarise(gene_list = list(unique(gene_symbol))) %>%
      ungroup() %>% select(gs_cat, gs_subcat, gs_name, gene_list)
    table$gene_list <- sapply(table$gene_list, toString)
    table_all <- rbind(table_all, table)
  }
}
write.csv(table_all, file = "gene_sets.csv", row.names = FALSE)