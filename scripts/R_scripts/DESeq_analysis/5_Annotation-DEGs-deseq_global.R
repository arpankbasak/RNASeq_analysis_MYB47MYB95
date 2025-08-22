rm(list = ls())

options(warn = 1,
       mc.cores = 6)

pkgs <- c("tidyverse", "GO.db", "parallel", "org.At.tair.db", "AnnotationDbi")
lapply(pkgs, require, character.only = TRUE)

# NOTE THAT RUNING ANALYSIS IS DONE IN EDGER

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
# load(paste(out.path, "/edger_objects.plotting.Rdata", sep = ""))
load(paste(out.path, "/LFC-km_cluster_edger-treatment.Rdata", sep = ""))

message("Conducting tidy annotation (BIOMART) ...")
logFC_P.annotated <- lfc_treatment
# logFC_P.annotated <- logFC_P

# Map all possible IDs
logFC_P.annotated$symbol = mapIds(org.At.tair.db,
                     keys=row.names(logFC_P.annotated), 
                     column="SYMBOL",
                     keytype="TAIR",
                     multiVals="list")
# logFC_P.annotated$enrtrez_id = mapIds(org.At.tair.db,
#                      keys=row.names(df), 
#                      column="ENTREZID",
#                      keytype="TAIR",
#                      multiVals="list")
logFC_P.annotated$gene_name =   mapIds(org.At.tair.db,
                     keys=row.names(logFC_P.annotated), 
                     column="GENENAME",
                     keytype="TAIR",
                     multiVals="list")
logFC_P.annotated$go_id =   mapIds(org.At.tair.db,
                     keys=row.names(logFC_P.annotated), 
                     column="GO",
                     keytype="TAIR",
                     multiVals="list")
logFC_P.annotated$enzyme =   mapIds(org.At.tair.db,
                     keys=row.names(logFC_P.annotated), 
                     column="ENZYME",
                     keytype="TAIR",
                     multiVals="list")
logFC_P.annotated$kegg_path =   mapIds(org.At.tair.db,
                     keys=row.names(logFC_P.annotated), 
                     column="PATH",
                     keytype="TAIR",
                     multiVals="list")
logFC_P.annotated$aracyc =   mapIds(org.At.tair.db,
                     keys=row.names(logFC_P.annotated), 
                     column="ARACYC",
                     keytype="TAIR",
                     multiVals="list")

# Tidy up the columns for annotations
logFC_P.annotated$gene_name <- sapply(logFC_P.annotated$gene_name, function(x) paste(x, collapse = ";")) %>% unlist(.) %>% as.character(.)
logFC_P.annotated$go_id <- sapply(logFC_P.annotated$go_id, function(x) paste(x, collapse = ";")) %>% unlist(.) %>% as.character(.)
logFC_P.annotated$enzyme <- sapply(logFC_P.annotated$enzyme, function(x) paste(x, collapse = ";")) %>% unlist(.) %>% as.character(.)
logFC_P.annotated$kegg_path <- sapply(logFC_P.annotated$kegg_path, function(x) paste(x, collapse = ";")) %>% unlist(.) %>% as.character(.)
logFC_P.annotated$aracyc <- sapply(logFC_P.annotated$aracyc, function(x) paste(x, collapse = ";")) %>% unlist(.) %>% as.character(.)
logFC_P.annotated$symbol <- sapply(logFC_P.annotated$symbol, function(x) paste(x, collapse = ";")) %>% unlist(.) %>% as.character(.)

message("DONE")

# Write output
save(list = c("logFC_P.annotated"),
     file = paste(out.path, "/gset_analysis/objects/edger-lfc_annotated.Rdata", sep = ""))
write.table(logFC_P.annotated, 
            paste(out.path, "/edger-lfc-p_annotated.txt", 
                  sep = ""), sep = "\t", quote = F)

# END OF SCRIPT
sessionInfo()