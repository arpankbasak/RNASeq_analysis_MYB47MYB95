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
clust.method <- "ward.D"

# Hierarchial clustering
idx <- str_detect(colnames(logFC_P), "^treatment") & str_detect(colnames(logFC_P), "_logFC$")
d <- 1 - cor(t(logFC_P[,idx]))
clust_gene_lfc <- hclust(as.dist(d), method = clust.method)
gene_clusters_lfc <- row.names(logFC_P[,idx])[clust_gene_lfc$order]

# Plot heatmap of clustering
clust_df <- as.data.frame(as.matrix(d)) %>% add_column(source = row.names(.), .before =1) %>%
    gather(key = "target", value = "vals", convert = FALSE, -source) %>%
    mutate(source = factor(source, levels = gene_clusters_lfc), 
           target = factor(target, levels = gene_clusters_lfc), 
           vals = saturate(1 - as.numeric(vals)))

clust_df %>%
    ggplot(aes(x = source, y=target)) +
    geom_raster(alpha = 1, aes(fill = vals)) +
    scale_fill_gradient2(low=gradient.low, 
                         mid = gradient.mid, 
                         high = gradient.high, 
                         na.value = gradient.na) +
    theme_AKB +
    theme(axis.text.x=element_blank(), 
          #strip.text.y=element_text(size = 10, face = "bold", angle = 180, hjust = 1, vjust = 0.5), 
          #strip.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_blank(),
          legend.text = element_text(size = 10)
         ) +
    labs(y="", x="", fill="") +
    ggsave(filename = paste(fig.path, "/COExp/LFC-sig[", alpha,"]-edger_treatment_COEXP_lfc-clustered[", clust.method,"]_heatmap.png", sep = ""),
           bg="transparent", width=4, height=3, limitsize=T, device = "png", dpi = 300)


# Kmeans clustering
n_centers=100
iter = 1000
#mds_obj <- cmdscale(d, k = n_centers, eig = T)
#eigens <- 1 - (mds_obj$eig/sum(mds_obj$eig))
aic_bic_list <- mclapply(1:n_centers, function(x){

    # x <- 4
    km_obj = kmeans(logFC_P[,idx], centers=x, iter.max = iter)
    aic_bic = kmeansAIC(km_obj)
    #repZ = repZ(km_obj, logFC_P[,idx], x)
    return(aic_bic)
    #hcc <- cutree(clust_gene_lfc, k = x)
    #temp <- silhouette(hcc, d)
    #spc <- specc(logFC_P[,idx], x)
}
, mc.cores = 16)

names(aic_bic_list) <- paste0("k_", as.character(1:n_centers))
aic_bic_df <- do.call(rbind, aic_bic_list)
conv <- aic_bic_df[,2]/sum(aic_bic_df[,2])
kx = 10
km_obj <- kmeans(logFC_P[,idx], centers= kx, iter.max = iter)
temp <- repZ(km_obj, logFC_P[,idx], kx)
clust_df$cluster_source <- as.factor(km_obj$cluster[match(clust_df$source, names(km_obj$cluster))])
clust_df$cluster_target <- as.factor(km_obj$cluster[match(clust_df$target, names(km_obj$cluster))])
# Plot Choice of clusters and stats for clustering --elbow plot

clust.method <- "km"

# Facet by clusters
clust_df %>%
    ggplot(aes(x = source, y=target)) +
    geom_raster(alpha = 1, aes(fill = vals)) +
    facet_grid(cluster_source ~ cluster_target, space = "free", scales = "free") +
    scale_fill_gradient2(low=gradient.low, 
                         mid = gradient.mid, 
                         high = gradient.high, 
                         na.value = gradient.na) +
    theme_AKB +
    theme(panel.spacing = unit(0.1, "lines"),
          axis.text.x=element_blank(), 
          strip.text.y = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5),
          strip.text.x = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.y=element_blank(),
          legend.text = element_text(size = 10)) +
    labs(y="", x="", fill="") +
    ggsave(filename = paste(fig.path, "/COExp/LFC-sig[", alpha,"]-edger_treatment_lfc-clustered_heatmap.png", sep = ""),
           bg="transparent", width=4, height=3, limitsize=T, device = "png", dpi = 300)

# Plot logFCP heatmap

           
# Consensus clustering
# idx <- str_detect(colnames(logFC_P), "^genotype") & str_detect(colnames(logFC_P), "_logFC$")
# d <- dist(logFC_P[,idx])
# clust_gene_lfc <- hclust(d, method = "ward.D")
# gene_clusters_lfc <- row.names(logFC_P[,idx])[clust_gene_lfc$order]

# Save clustering objects
save(list = c("clust_gene_lfc", "gene_clusters_lfc"),
     file = paste(out.path, "/LFC-cluster_edger-genotype.Rdata", sep = ""))


# END OF SCRIPT
sessionInfo()
