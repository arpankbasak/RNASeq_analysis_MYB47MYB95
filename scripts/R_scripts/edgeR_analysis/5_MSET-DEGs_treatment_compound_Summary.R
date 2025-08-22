# Script for MSET analysis using differential gene as input across genotype
# This script is written to fetch compound information for integrative analysis of metabolome and RNASeq data
# @Arpan Kumar Basak

# Cleanup
rm(list = ls())

options(warn = 1,
       mc.cores = 6)

pkgs <- c("tidyverse", "fgsea", "GO.db", "parallel", "KEGGREST")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
load(paste(out.path, "gset_analysis/objects/edeger-GO-Pathway_genotype.Rdata", sep = ""))

# Link compounds and pathway associaktions
compound_list <- KEGGREST::keggLink("pathway", "compound")

for(i in 1:nrow(path_df)){
    
    idx <- which(compound_list %in% path_df$kegg_path_id[i])
    path_df$compound[i] <- paste(unique(names(compound_list)[idx]), collapse = ";")
}

# Make empty data frame
compound_df <- data.frame(compound_id = unique(names(compound_list)), stringsAsFactors = FALSE) %>%
    add_column(exact_mass = NA, 
               formula = NA,
               entry = NA)
message("Fetching compounds for pathway that is differentially expressed ...")

# Fetch information of compounds as EXACT MASS
compound_df$exact_mass <- sapply(compound_df$compound_id, function(x) {
    KEGGREST::keggGet(x)[[1]]$EXACT_MASS
}, USE.NAMES = FALSE)

# Fetch formula
compound_df$formula <- sapply(compound_df$compound_id, function(x) {
    KEGGREST::keggGet(compound_df$compound_id[i])[[1]]$FORMULA
}, USE.NAMES = FALSE)

# Fetch enry information
compound_df$entry <- sapply(compound_df$compound_id, function(x) {
    unname(KEGGREST::keggGet(compound_df$compound_id[i])[[1]]$ENTRY)
}, USE.NAMES = FALSE)

# Save output
save(list = "compound_df", 
     file = paste(out.path, "MSET_edeger-GO-KEGG_Pathway_genotype.Rdata", sep = ""))

message("DONE")

# END OF SCRIPT
sessionInfo()