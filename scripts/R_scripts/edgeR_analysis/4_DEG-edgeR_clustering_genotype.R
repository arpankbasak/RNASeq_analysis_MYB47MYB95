# Scripts for clustering differential genes from EDGER output
# @Arpan Kumar Basak

# Cleanup
rm(list = ls())

options(warn = 1,
       mc.cores = 8)

pkgs <- c("tidyverse", "kernlab", "cluster", "parallel")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
load(paste(out.path, "/edger_objects.plotting.Rdata", sep = ""))

# Cluster genes on the basis of LFC- expression

# Hierarchial clustering
idx <- str_detect(colnames(logFC_P), "^genotype") & str_detect(colnames(logFC_P), "_logFC$")
d <- dist(logFC_P[,idx])
clust_gene_lfc <- hclust(d, method = "ward.D")
gene_clusters_lfc <- row.names(logFC_P[,idx])[clust_gene_lfc$order]

# Kmeans clustering
n_centers=100
iter = 100
#mds_obj <- cmdscale(d, k = n_centers, eig = T)
#eigens <- 1 - (mds_obj$eig/sum(mds_obj$eig))
aic_bic_list <- mclapply(1:n_centers, function(x){

    # x <- 4
    km_obj <- kmeans(logFC_P[,idx], centers=x, iter.max = iter)
    aic_bic <- kmeansAIC(km_obj)
    repZ <- repZ(km_obj, logFC_P[,idx], x)
    #hcc <- cutree(clust_gene_lfc, k = x)
    #temp <- silhouette(hcc, d)
    #spc <- specc(logFC_P[,idx], x)
 
  }

, mc.cores = 12)

names(aic_bic_list) <- paste0("k_", as.character(1:n_centers))
aic_bic_df <- do.call(rbind, aic_bic_list)
km_obj <- kmeans(logFC_P[,idx], centers=86, iter.max = iter)
temp <- repZ(km_obj, logFC_P[,idx], 86)

# Consen clustering
idx <- str_detect(colnames(logFC_P), "^genotype") & str_detect(colnames(logFC_P), "_logFC$")
d <- dist(logFC_P[,idx])
clust_gene_lfc <- hclust(d, method = "ward.D")
gene_clusters_lfc <- row.names(logFC_P[,idx])[clust_gene_lfc$order]

# Save clustering objects
save(list = c("clust_gene_lfc", "gene_clusters_lfc"),
     file = paste(out.path, "/LFC-cluster_edger-genotype.Rdata", sep = ""))


# END OF SCRIPT
sessionInfo()
