rm(list = ls())

options(warn = 1,
       mc.cores = 6)

pkgs <- c("tidyverse", "GO.db", "parallel", "org.At.tair.db", "AnnotationDbi")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
load(paste(out.path, "/DESeq2_objects.Rdata", sep = ""))

message("Conducting tidy annotation (BIOMART) ...")
df <- logFC_P

# Map all possible IDs
df$symbol = mapIds(org.At.tair.db,
                     keys=row.names(df), 
                     column="SYMBOL",
                     keytype="TAIR",
                     multiVals="list")
df$gene_name = mapIds(org.At.tair.db,
                     keys=row.names(df), 
                     column="ENTREZID",
                     keytype="TAIR",
                     multiVals="list")
df$name =   mapIds(org.At.tair.db,
                     keys=row.names(df), 
                     column="GENENAME",
                     keytype="TAIR",
                     multiVals="list")
df$name =   mapIds(org.At.tair.db,
                     keys=row.names(df), 
                     column="GO",
                     keytype="TAIR",
                     multiVals="list")
df$name =   mapIds(org.At.tair.db,
                     keys=row.names(df), 
                     column="ENZYME",
                     keytype="TAIR",
                     multiVals="list")
df$name =   mapIds(org.At.tair.db,
                     keys=row.names(df), 
                     column="PATH",
                     keytype="TAIR",
                     multiVals="list")

# Gene information
df$gene_name <- sapply(row.names(df), 
                       function(x) paste(unique(call_biomart(x,
                                   columns = c("external_gene_name"))$external_gene_name), 
                      collapse = ";"))
df$chromosomhead()e <- sapply(row.names(df), 
                        function(x) paste(unique(call_biomart(x,
                                   columns = c("chromosome_name"))$chromosome_name), 
                                          collapse = ""))
df$start_position <- sapply(row.names(df), 
                            function(x) paste(unique(call_biomart(x,
                                   columns = c("start_position"))$start_position), 
                                              collapse = ""))
df$end_position <- sapply(row.names(df), 
                          function(x) paste(unique(call_biomart(x,
                                   columns = c("end_position"))$end_position), 
                                            collapse = ""))
df$strand <- sapply(row.names(df), 
                    function(x) paste(unique(call_biomart(x,
                                   columns = c("strand"))$strand), 
                                      collapse = ""))


# GO information
df$go_id <- sapply(row.names(df), function(x)
    paste(unique(call_biomart(x, columns = c("go_id"))$go_id), 
          collapse = ";"))

df$go_definition <- sapply(row.names(df), function(x)
    paste(unique(call_biomart(x, columns = c("definition_1006"))$definition_1006), 
          collapse = ";"))

# Plant pathway information
df$po_ids <- sapply(row.names(df), function(x)
    paste(call_biomart(x,columns = c("po_id"))$po_id, 
          collapse = ";"))

df$po_definition <- sapply(row.names(df), function(x)
    paste(unique(call_biomart(x,
                       columns = c("po_definition_1006"))$po_definition_1006), 
          collapse = ";"))
                           
df <- cbind.data.frame(df, logFC_P.cook[row.names(df), c("dispOutlier", "maxCooks")])

logFC_P.annotated <- df

message("DONE")

# Write output
save(list = c("logFC_P.annotated"),
     file = paste(out.path, "/gset_analysis/objects/edger-lfc_annotated.Rdata", sep = ""))
write.table(logFC_P.annotated, 
            paste(out.path, "/edger-lfc-p_annotated.txt", sep = ""), sep = "\t")

# END OF SCRIPT
sessionInfo()