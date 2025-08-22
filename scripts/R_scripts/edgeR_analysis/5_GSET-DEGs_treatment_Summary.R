rm(list = ls())

# TIDY Summary of up and down regulated genes in GO and Pathway terms
options(warn = 1,
       mc.cores = 6)

pkgs <- c("tidyverse", "limma", "GO.db", "parallel", "KEGGREST")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
load(paste(out.path, "gset_analysis/objects/edeger-biomart_treatment_up.Rdata", sep = ""))

df_up <- data.frame()

message("Fetching Gene Ontology terms for upregulated genes...")
for(i in 1:length(bm_res_upregulated)){
    df_up <- rbind(df_up, distinct(bm_res_upregulated[[i]], go_id, .keep_all = TRUE))
}

# Tidy for GO representation for conditin specific
go_df_up <- df_up %>% dplyr::select(ensembl_gene_id, external_gene_name, go_id) %>%
        group_by(go_id) %>% 
        summarise(gene_id = paste(unique(ensembl_gene_id), collapse = ";"),
                  gene_name = paste(unique(external_gene_name), collapse = ";")) %>%
        mutate(regulation = "up") %>%
        data.frame(., stringsAsFactors = FALSE)

# Call GO for annotation store the symbol
path_df_up <- call_GO(go_df_up$go_id, by = "GO") %>% 
            mutate(kegg_path_id = ifelse(!is.na(PATH), 
                                         paste("path:map", PATH, sep = ""), "unknown")) %>%
            distinct(kegg_path_id, .keep_all = TRUE) %>%
            mutate(regulation = "up") %>%
            data.frame(., stringsAsFactors = FALSE) %>% na.omit(.)


# Save GO objects and Path objetcs as well as data fremes
rm(list = c("bm_res_upregulated"))
df_up$regulation <- "up"

message("Fetching Gene Ontology terms for downregulated genes...")
load(paste(out.path, "gset_analysis/objects/edeger-biomart_treatment_down.Rdata", sep = ""))

df_down <- data.frame()
for(i in 1:length(bm_res_downregulated)){
    df_down <- rbind(df_down, distinct(bm_res_downregulated[[i]], go_id, .keep_all = TRUE))
}

# Tidy for GO representation
go_df_down <- df_down %>% dplyr::select(ensembl_gene_id, external_gene_name, go_id) %>%
        group_by(go_id) %>% 
        summarise(gene_id = paste(unique(ensembl_gene_id), collapse = ";"),
                  gene_name = paste(unique(external_gene_name), collapse = ";")) %>%
        mutate(regulation = "down") %>%
        data.frame(., stringsAsFactors = FALSE)

# Call GO for annotation store the symbol
path_df_down <- call_GO(go_df_down$go_id, by = "GO") %>% 
            mutate(kegg_path_id = ifelse(!is.na(PATH), 
                                         paste("path:map", PATH, sep = ""), "unknown")) %>%
            distinct(kegg_path_id, .keep_all = TRUE) %>%
            mutate(regulation = "down") %>%
            data.frame(., stringsAsFactors = FALSE) %>% na.omit(.)

go_df <- rbind(go_df_up, go_df_down)
path_df <- rbind(path_df_up, path_df_down)

# Save GO objects and Path objetcs as well as data fremes
save(list = c("go_df", "path_df"), 
     file = paste(out.path, "gset_analysis/objects/edeger-GO-Pathway_treatment.Rdata", sep = ""))
rm(list = c("bm_res_downregulated"))
df_down$regulation <- "down"

biomart_summary_table <- rbind(df_up, df_down)
save(list = c("biomart_summary_table"), 
     file = paste(out.path, "gset_analysis/objects/edeger-biomart-summary_treatment.Rdata", sep = ""))
write.table(biomart_summary_table, 
            file = paste(out.path, "gset_analysis/objects/edeger-biomart-summary_treatment.txt", sep = ""),
            sep = "\t")

# ocncatenate both up and down regulated once

message("DONE")

# END OF SCRIPT
sessionInfo()