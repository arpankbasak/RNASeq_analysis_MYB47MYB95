# Script for doing exploratory analysis
# @Arpan Kumar Basak
# Cleanup
rm(list = ls())

options(warn = 1,
       mc.cores = 8)

pkgs <- c("tidyverse", "biocparallel", "DESeq2")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
#edata <- read.table(paste(out.path, "/summarised_exp/edata_rs.txt", sep = ""),
#                       header = TRUE, as.is = TRUE)
load(paste(out.path, "/summarised_exp/salmon_summarized_experiments_genelevel.Rdata", sep = ""))
load(paste(out.path, "edger_objects.plotting.Rdata", sep = ""))
load(paste(out.path, "/DESeq2_objects.Rdata", sep = ""))
idx <- row.names(log2cpm)
# loop to plot MDS for all matrix
list.norm <- list(RLD = assay(rld)[idx,], 
                  CPM = log2cpm, 
                  beta = apply(tximport_sum$counts[idx,], 2, function(x) x/sum(x)))

mclapply(1:length(list.norm), function(i) {
    
    # i = 2
    x <- list.norm[[i]]
    id <- names(list.norm)[i]
    
    d <- dist(t(x))
    clust_sample <- hclust(d, method = "ward.D")
    sample_clusters <- colnames(x)[clust_sample$order]
    mds_obj <- cmdscale(d, k = 3, eig = TRUE)
    eig <- mds_obj$eig/sum(mds_obj$eig)
    idx <- match(design$sample_name, row.names(mds_obj$points))
    df <- cbind(design, 
            mds1 = mds_obj$points[design$sample_name,1], 
            mds2 = mds_obj$points[design$sample_name,2], 
            mds3 = mds_obj$points[design$sample_name,3])

    message(paste0("Spilling MDS Plot for ", id," ..."))
    df %>% 
    # mutate(sample_name = factor(sample_name, 
    #                                    levels = sample_clusters)) %>%
            ggplot(aes(x = saturate(mds1), y = saturate(mds2))) +
            ggtitle(paste0(id)) +
            # geom_point(alpha = 0.8, aes(shape = treatment, colour = genotype)) +
            geom_text(aes(colour = genotype, label = sample_name), 
                      # nudge_x = 0.5, nudge_y = 0.5,
                      check_overlap = FALSE, 
                      show.legend = FALSE, 
                      size = 1) +
            scale_colour_manual(values = genotype$col.idx, labels = label.idx) +
            scale_shape_manual(values = treatment$pch.idx) +
            stat_ellipse(level = 0.95, type = "t", aes(group = treatment, linetype = treatment)) +
            scale_linetype_manual(values = treatment$lty.idx) +
            geom_hline(yintercept = 0.0, linetype = 4, colour = "lightgrey") +
            geom_vline(xintercept = 0.0, linetype = 4, colour = "lightgrey") +
            # xlim(c(-1,1)) +
            # ylim(c(-1,1)) +
            theme_AKB +
            labs(x = paste("MDS 1: ", round(eig[1], 2) * 100, " % Var. Exp.", sep = ""), 
                 y = paste("MDS 2: ", round(eig[2], 2) * 100, " % Var. Exp.", sep = ""), 
                 colour = "",
                 shape = "", linetype = "") +
            theme(legend.title = element_blank(),
                  legend.text = element_text(size = 6)) +
            ggsave(paste(fig.path, "/EDA_DEGs_MDS12_", id,"_Plot.png", sep = ""), 
                   dpi = 300, device = "png", 
                   height = 5, units = "in", 
                   width = 6, limitsize = F)
    
     df %>% 
     # mutate(sample_name = factor(sample_name, levels = sample_clusters)) %>%
            ggplot(aes(x = saturate(mds2), y = saturate(mds3))) +
            ggtitle(paste0(id)) +
            # geom_point(alpha = 0.8, aes(shape = treatment, colour = genotype)) +
            geom_text(aes(colour = genotype, label = sample_name), 
                      # nudge_x = 0.5, nudge_y = 0.5,
                      check_overlap = FALSE, 
                      show.legend = FALSE, 
                      size = 1) +
            scale_colour_manual(values = genotype$col.idx, labels = label.idx) +
            scale_shape_manual(values = treatment$pch.idx) +
            stat_ellipse(level = 0.95, type = "t", aes(group = treatment, linetype = treatment)) +
            scale_linetype_manual(values = treatment$lty.idx) +
            geom_hline(yintercept = 0.0, linetype = 4, colour = "lightgrey") +
            geom_vline(xintercept = 0.0, linetype = 4, colour = "lightgrey") +
            # xlim(c(-1,1)) +
            # ylim(c(-1,1)) +
            theme_AKB +
            labs(x = paste("MDS 2: ", round(eig[2], 2) * 100, " % Var. Exp.", sep = ""), 
                 y = paste("MDS 3: ", round(eig[3], 2) * 100, " % Var. Exp.", sep = ""), 
                 colour = "",
                 shape = "", linetype = "") +
            theme(legend.title = element_blank(),
                  legend.text = element_text(size = 6)) +
            ggsave(paste(fig.path, "/EDA_DEGs_MDS23_", id,"_Plot.png", sep = ""), 
                   dpi = 300, device = "png", 
                   height = 5, units = "in",
                   width = 6, limitsize = F)
        
        k <- 10
        d <- dist(x)
        clust_genes <- hclust(d, method = "ward.D")
        gene_clusters <- colnames(x)[clust_sample$order]
        
        # Add clusters by kmeans, agglomerative clustering and spectral clustering
#         hc <- cutree(clust_genes, k = k)
#         km <- kmeans(x, centers = k)
#         spc <- specc()
    
        df <- as.data.frame(x) %>% 
            add_column(gene_ids = row.names(x), .before = 1) %>%
            gather(key = "sample_name", value = "counts", convert = FALSE, -gene_ids) %>%
            mutate(counts = ifelse(counts != 0, log2(counts + 1), log2(1)),
                  sample_name = factor(sample_name, levels = sample_clusters),
                  gene_ids = factor(gene_ids, levels = gene_clusters)) # Add a pseudocounts

        df$genotype <- design$genotype[match(df$sample_name, design$sample_name)]
        df$treatment <- design$treatment[match(df$sample_name, design$sample_name)]

        df$genotype <- factor(df$genotype, levels = genotype$name)
        df$treatment <- factor(df$treatment, levels = treatment$name)

        levels(df$genotype) <- genotype$short

        # Heatmap for counts
        # df %>% ggplot(aes(x = sample_name, 
        #                   y = gene_ids, 
        #                   fill = saturate(counts))) +
        #     ggtitle(paste0(id)) +
        #     geom_raster(alpha = 1) +
        #     facet_grid(~ treatment + genotype, scales = "free",
        #                switch = "x", 
        #                labeller = labeller(label.idx)) +
        #     scale_fill_gradient(na.value = "darkgrey", 
        #                         low = "black",
        #                        high = "lightyellow") +
        #     theme_AKB +
        #     theme(axis.text.x = element_blank(),
        #           axis.text.y = element_blank(),
        #           strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        #           legend.text = element_text(size = 6)) +
        #     labs(x = "", y = "",
        #         fill = "log2[counts]") +
        #     ggsave(paste(fig.path, "/EDA_heatmap_",id,"_genes.png", sep = ""), 
        #            dpi = 300, 
        #            device = "png", 
        #            units = "in",
        #            height = 8, 
        #            width = 4, limitsize = F)

        
    
}, mc.cores = 4)
                               

message("DONE")

# END OF SCRIPT
sessionInfo()
