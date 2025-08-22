# Script for GSET analysis using MROAST and CAMERA pipelines
# @Arpan Kumar Basak

# Cleanup
rm(list = ls())
# NOTE::
# Gene set here is the list of gese that were significantly detected across genotype within treatment conditions
options(warn = 1,
       mc.cores = 8)

pkgs <- c("tidyverse", "biocparallel", "DESeq2", "limma", "edgeR")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
edata <- read.table(paste(out.path, "/summarised_exp/edata.filtered_rs.txt", sep = ""),
                       header = TRUE, as.is = TRUE)
load(paste(out.path, "/DESeq2_objects.Rdata", sep = ""))
design <- design %>% mutate(random = paste("B", random, sep = ""), 
                            seq_run = paste("B", seq_run, sep = ""))

row.names(design) <- design$sample_name
load(paste(out.path, "gset_analysis/objects/edeger-biomart-summary_genotype.Rdata", sep = ""))

# Build FDATA -- on the basis of up and down regulation of genes
fdata <- biomart_summary_table %>% 
            dplyr::select(go_id, ensembl_gene_id, external_gene_name, regulation) %>% 
            group_by(regulation, go_id) %>% 
            summarise(gene_ids = paste(unique(ensembl_gene_id), collapse = ";"), 
                      gene_name = paste(unique(external_gene_name), collapse = ";")) %>%
            data.frame(., stringsAsFactors = FALSE)

# BUILD Annotation data list for MROAST and CAMERA
adata <- list()
for(i in unique(fdata$go_id)){
    adata[[i]] <- unlist(str_split(fdata$gene_id[fdata$go_id == i], ";"), use.names = FALSE)
    adata[[i]] <- adata[[i]][!is.na(adata[[i]])]
}

# Using ARACYC database
aracyc <- call_GO(unique(fdata$go_id), "ARACYC", by = "GO") %>% 
            group_by(GO) %>% 
            summarise(aracyc = paste(unique(ARACYC), collapse = ";")) %>% 
            dplyr::select(GO, aracyc) %>%
            mutate(aracyc = str_replace_all(aracyc, "NA", "unknown")) %>%
            data.frame(.)

# Load GLM objects - Using the same contrast as used in GLM analysis at gene level
load(paste(out.path, "edger_objects.GLM.Rdata", sep = ""))

# Build indeces for the analysis
idx <- ids2indices(adata, row.names(de))
idy <- unlist(lapply(idx, length))
idy <- which(between(idy, 3, 200)) # Specify range of number of genes to be analysed a.k.a set Target genes

# Conduct MROAST --Targetted
message("Conducting geneset analysis using GO-tags and MROAST ...")

mroast.list <- lapply(contrast.names, function(x) 
    mroast(de, idx[idy], model, contrast=contrast.mat[, which(contrast.names == x)], nrot = 1000))
# make tidy MROAST
mroast.GO.list <- lapply(1:n, function(x) {
  table <- as.data.frame(mroast.list[[x]][,c("NGenes", "Direction", "PValue", "FDR")])
  table <- table %>% add_column(go_id = row.names(table), .before = 1) %>%
    #dplyr::select(go_id, NGenes, PValue, FDR) %>% 
    dplyr::mutate(contrast = contrast.names[x])
  #colnames(table)[-1] <- paste(contrast.names[x], colnames(table)[-1], sep="_")
  return(table)
})

# Tidy up the table
mroast.GO.long <- do.call(rbind, mroast.GO.list) 
idx <- match(mroast.GO.long$go_id, fdata$go_id)
mroast.GO.long$gene_name <- fdata$gene_name[idx]
idx <- match(mroast.GO.long$go_id, aracyc$GO)
mroast.GO.long$aracyc <- aracyc$aracyc[idx]
write.table(mroast.GO.long, 
            paste(out.path, "mroast-edger_GSET_long_genotype.txt", sep = ""), sep = "\t")

message("DONE")                      

# Conduct CAMERA analysis --Untargetted
message("Conducting geneset analysis using GO-tags and CAMERA ...")
                      
idx <- ids2indices(adata, row.names(de))
# Conduct CAMERA
camera.list <- lapply(contrast.names, function(x) 
    camera(de, idx, model, contrast=contrast.mat[, which(contrast.names == x)]))
# make tidy CAMERA
camera.GO.list <- lapply(1:n, function(x) {
  table <- as.data.frame(camera.list[[x]])
  table <- table %>% add_column(go_id = row.names(table), .before = 1) %>%
    dplyr::mutate(contrast = contrast.names[x])
  return(table)
})

camera.GO.long <- do.call(rbind, camera.GO.list) 
idx <- match(camera.GO.long$go_id, fdata$go_id)
camera.GO.long$gene_name <- fdata$gene_name[idx]
idx <- match(camera.GO.long$go_id, aracyc$GO)
camera.GO.long$aracyc <- aracyc$aracyc[idx]

write.table(camera.GO.long, 
            paste(out.path, "camera-edger_GSET_long_genotype.txt", sep = ""), sep = "\t")
                      
save(list = c("mroast.GO.long", "camera.GO.long"), 
     file = paste(out.path, "GSET-GO_edger_objects-genotype.plotting.Rdata", sep = ""))
                      
message("DONE")

# END OF SCRIPT
sessionInfo()