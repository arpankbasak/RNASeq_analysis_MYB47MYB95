rm(list = ls())

# Normalised by WT treatment expression
options(warn = 1,
       mc.cores = 8)

pkgs <- c("tidyverse", "VennDiagram", "DESeq2", "grid")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)
set.seed(seed)

# Raad in data
load(paste(out.path, "/edger_objects.plotting.Rdata", sep = ""))
# load(paste(out.path, "/LFC-km_cluster_edger-treatment.Rdata", sep = ""))
load(paste(out.path, "/gset_analysis/objects/edger-lfc_annotated.Rdata", sep = ""))

# Untargetted effect of treatment within genotype using k-means clustering
# idx <- str_detect(colnames(logFC_P), "^genotype_")
# logFC_P <- logFC_P[,idx]
# logFC_P <- cbind.data.frame(logFC_P[,str_detect(colnames(logFC_P), "_logFC$")] - logFC_P[,"treatment_WT_JA_logFC"], logFC_P[,!str_detect(colnames(logFC_P), "_logFC$")])

# lfc_df <- logFC_P %>%
#     add_column(geneids = row.names(.), .before = 1) %>%
#     # cbind.data.frame(., cluster_k = lfc_treatment$cluster_k[match(.$geneids, row.names(lfc_treatment))]) %>%
#     gather(key = "interaction", value = "vals", convert = FALSE, -geneids) %>%
#     separate(interaction, into = c("x", "x1", "genotype", "treatment", "x2"), 
#              convert = F, sep = "_") %>%
#     dplyr::select(-x, -x1) %>%
#     spread(key = "x2", value = "vals", convert = FALSE) %>%
#     mutate(PValue = replace_na(PValue, 1)) %>%
#     mutate(genotype = factor(genotype, 
#                              levels = c("WT", "myb", "myc.tKO")),
#            sig = as.factor(ifelse(PValue < 0.05, "Significant", "NS")),
#            # treatment = as.factor(treatment$name[match(treatment, design$sample_name)]),
#            geneids = as.factor(geneids)
#           )

lfc_annotated <- cbind.data.frame(logFC_P, logFC_P.annotated)
# row.names(lfc_annotated) <- row.names(logFC_P)
write.table(lfc_annotated, "./statistics/annotated_genotype_comparison_normalised.txt", sep = "\t", row.names = TRUE,  quote = FALSE)

# Plot of tracer
# logFC_P[target,] %>% 
#     add_column(geneids = row.names(.), .before = 1) %>%
#     gather(key = "interaction", value = "vals", convert = FALSE, -geneids) %>%
#     separate(interaction, into = c("x", "genotype", "x1", "x2"), 
#              convert = F, sep = "_") %>%
#     spread(key = "x2", value = "vals", convert = FALSE) %>%
#     mutate(PValue = replace_na(PValue, 1)) %>%
#     mutate(genotype = factor(genotype, 
#                              levels = c("WT", "myb", "myc.tKO")),
#            sig = as.factor(ifelse(PValue < 0.05, "Significant", "NS")),
# #            treatment = as.factor(treatment$name[match(treatment, design$sample_name)]),
#            geneids = as.factor(geneids)
#           ) %>%
#     ggplot(aes(y = logFC, x = genotype)) +
#     ggtitle(paste("logFC vs DW", sep = "")) +
#     geom_bar(fill = "transparent", 
#              aes(colour = sig),
#              stat = "identity"
#             ) +
#     facet_wrap(geneids ~ .
# #                labeller = labeller(label.idx)
#               ) +
#     scale_colour_manual(values = c(`Significant` = "black", 
#                                  `NS` = "grey")
#                      ) +
# #     scale_shape_manual(values = genotype$pch.idx, labels = label.idx) +
#     theme_AKB +
#     theme(axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5, size = 10),
#           axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
#           strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
#           legend.text = element_text(size = 6)
#          ) +
#     labs(x = "", colour = "", shape = "", y = "log2FC[counts]") +
#     ggsave(paste(fig.path, "/lfc_diagnosis_hallmark.png", sep = ""),
#            dpi = 600,
#            device = "png", 
#            bg = "transparent",
#            height = 7,
#            width = 14, limitsize = F)

# Remove genes that donot pass Procedure for MHT
# before_MHT <- nrow(logFC_P)
# logFC_P <- logFC_P %>% na.omit(.)
# after_MHT <- nrow(logFC_P)

# Cluster genes on the basis of LFC- expression
if(file.exists(paste(out.path, "/LFC-cluster_edger-treatment.Rdata", sep = ""))){

  load(paste(out.path, "/LFC-cluster_edger-treatment.Rdata", sep = ""))

}else {
  
  idx <- str_detect(colnames(logFC_P), "^treatment") & str_detect(colnames(logFC_P), "_logFC$")
  d <- 1 - cor(t(logFC_P[,idx]))
  clust_gene_lfc <- hclust(as.dist(d), method = "ward.D")
  gene_clusters_lfc <- row.names(logFC_P[,idx])[clust_gene_lfc$order]
  save(list = c("clust_gene_lfc", "gene_clusters_lfc"),
       file = paste(out.path, "/LFC-cluster_deseq-treatment.Rdata", sep = ""))
}



# Plot cluster representives of the gene
# feature_dendrogram(clust_gene_lfc) +
#   ggsave(filename = paste0(fig.path, "/DEGs/DEFs_deseq_hclusters_(p_", alpha,"_permil_LFC: ",threshold,")_genotype_dendrogram_lfc_sorted.png"),
#     bg="transparent", width=3.5, height=14, limitsize=F, device = "png", dpi = 300)


# Differential Gene Data carpenting
idx <- str_detect(colnames(logFC_P), "^treatment")
df_lfc <- logFC_P[,idx] %>% 
            add_column(gene_id = row.names(logFC_P[,idx]), .before = 1) %>%
            gather(key = "key", value = "val", convert = FALSE, -gene_id) %>%
            separate(key, into = c("x1", "genotype", "x2", "col"), 
                     convert = FALSE, sep = "_") %>%
            dplyr::select(-x1,-x2) %>%
            spread(key = "col", value = "val", convert = FALSE) %>%
            mutate(PValue = replace_na(PValue, 1)) %>%
            mutate(gene_id = factor(gene_id, levels = gene_clusters_lfc)) %>%
            data.frame(.)

# Make factors
df_lfc$genotype <- factor(df_lfc$genotype, levels= genotype$short)
df_lfc$sig <- as.factor(ifelse(df_lfc$PValue < alpha & abs(df_lfc$logFC) > 0.5, 1 * sign(df_lfc$logFC), 0))


message("Spilling LFC heatmap ...")
# # PLot Canvas for LFC
# df_lfc %>% ggplot(aes(x = genotype, y = gene_id, fill = saturate(logFC))) +
#     geom_raster(alpha = 1) +
#     scale_fill_gradient2(low = gradient.low, 
#                          mid = gradient.mid,
#                          high = gradient.high, 
#                          na.value = gradient.na) +
#     theme_AKB +
#     theme(panel.spacing = unit(0.2, "lines"),,
#       axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
#           axis.text.y = element_blank()) +
#     labs(x = "",
#          y = "",
#          fill = "log2(Fold_Change)") +
#     ggsave(paste(fig.path, "/DEGs/LFC-deseq_treatment_lfc-clustered_heatmap.png", sep = ""), 
#            dpi = 600, 
#            device = "png", 
#            units = "in",
#            height = 10, 
#            width = 4, 
#            limitsize = F, bg = "transparent")

# message("Spilling LFC heatmap gold standards ...")
# # Those which are significant
# df_lfc %>% ggplot(aes(x = genotype, y = gene_id, fill = sig)) +
#     geom_raster(alpha = 1) +
#     scale_fill_manual(values = c(`0` = gradient.na,
#                                  `1` = gradient.sig.high, 
#                                  `-1` = gradient.sig.low),
#                      labels = c(`0` = "Not Significant",
#                                  `1` = "Upregulated", 
#                                  `-1` = "Downregulated")) +
#     theme_AKB +
#     theme(panel.spacing = unit(0.2, "lines"),,
#       axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
#           axis.text.y = element_blank()) +
#     labs(x = "",
#          y = "",
#          fill = paste("FDR < ", alpha, sep = "")) +
#     ggsave(paste(fig.path, "/DEGs/LFC-sig[", alpha,"]-deseq_treatment_lfc-clustered_heatmap.png", sep = ""), 
#            dpi = 600, 
#            device = "png", 
#            height = 21, 
#            width = 7, limitsize = F, bg = "transparent")

load(paste(out.path, "/log.normalised.matrix.Rdata", sep = ""))
df_lfc$high_sig <- ifelse(-log10(df_lfc$PValue) > 15, df_lfc$gene_id, NA)
idx <- match(df_lfc$gene_id, row.names(logFC_P.annotated))
df_lfc <- cbind.data.frame(df_lfc, 
  logFC_P.annotated[idx, !str_detect(colnames(logFC_P.annotated), "^treatment")])

cluster_mat <- lfc_treatment %>%
gather(key = "key", value = "vals", -cluster_k) %>%
group_by(cluster_k, key) %>%
summarise(vals = mean(vals)) %>%
spread(key = key, value = vals, fill = 0) %>%
data.frame(.)

d <- 1-cor(t(cluster_mat[,-1]))
hc <- hclust(as.dist(d))
cluster_sorted <- cluster_mat[,1][hc$order]

# Clustering of genes for Log expression and 

# Heatmap for counts
# df %>% mutate(gene_ids = factor(as.character(gene_ids), levels = gene_clusters_lfc)) %>%
#     ggplot(aes(x = sample_name, 
#                   y = gene_ids, 
#                   fill = saturate(counts))) +
#     geom_raster(alpha = 1) +
#     facet_grid(~ genotype + treatment,
#                switch = "x", space = "free", scale = "free",
#                labeller = labeller(label.idx)) +
#     scale_fill_gradient(na.value = "darkgrey", 
#                         low = "black", 
#                        high = "yellow") +
#     theme_AKB +
#     theme(panel.spacing = unit(0.2, "lines"),
#       axis.text.x = element_blank(),
#           axis.text.y = element_blank(),
#           strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
#           legend.text = element_text(size = 6)) +
#     labs(x = "", y = "",
#         fill = "log2[counts]") +
#     ggsave(paste(fig.path, "/EDA_heatmap_treatment_lfc_clustered.png", sep = ""), 
#            dpi = 600, 
#            device = "png", 
#            height = 21, 
#            width = 7, limitsize = F, bg = "transparent")



# Plot of untargetted K-Means clustering for the effect of JA treatment within genotype
# strip for kmeans
n_clusters <- lfc_treatment %>%
group_by(cluster_k) %>%
summarise(
  n_genes = n(), 
  up_wt = sum(sign(treatment_WT_JA_logFC) == 1),
  down_wt = sum(sign(treatment_WT_JA_logFC) == -1),
  up_myb = sum(sign(treatment_myb_JA_logFC) == 1),
  down_myb = sum(sign(treatment_myb_JA_logFC) == -1),
  up_myc = sum(sign(treatment_myc.tKO_JA_logFC) == 1),
  down_myc = sum(sign(treatment_myc.tKO_JA_logFC) == -1)) %>%
data.frame(., stringsAsFactors = FALSE)


(strip_cluster <- df_lfc %>%
mutate(
  gene_id = factor(gene_id, levels = gene_clusters_lfc),
  cluster_k = factor(cluster_k, levels = cluster_sorted),
  cluster_dens = log10(n_clusters$n_genes[match(as.character(.$cluster_k), n_clusters$cluster_k)])
  ) %>%
    ggplot(aes(x = "", y = gene_id, 
                  fill = (cluster_dens))) +
    geom_raster(alpha = 1) +
    # geom_tile(fill = NA, width = 0.95, height = 0.95, aes(colour = sig)) +
    facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
    scale_fill_gradient(low = "lightgrey", 
                         high = "darkgrey", 
                         na.value = gradient.na) +
    # scale_colour_manual(values = c(`NS` = NA, `Significant` = "black")) +
    theme_AKB +
    theme(panel.spacing = unit(0.2, "lines"),
          strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
          axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
          axis.text.y = element_blank()) +
    labs(x = "",
         y = "",
         fill = "")) +
    ggsave(paste(fig.path, "/DEGs/km_clusters_heatmap.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 20, 
           width = 4, limitsize = F, bg = "transparent")

(lfc_hmap <- df_lfc %>%
 mutate(gene_id = factor(gene_id, levels = gene_clusters_lfc),
        cluster_k = factor(cluster_k, levels = cluster_sorted)) %>%
     ggplot(aes(x = genotype, y = gene_id, 
                   fill = saturate(logFC))) +
     geom_raster(alpha = 1) +
     # geom_tile(fill = NA, width = 0.95, height = 0.95, aes(colour = sig)) +
     facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
     scale_fill_gradient2(low = gradient.low, 
                          mid = gradient.mid,
                          high = gradient.high, 
                          na.value = gradient.na) +
     # scale_colour_manual(values = c(`NS` = NA, `Significant` = "black")) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
           strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
           axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
           axis.text.y = element_blank()) +
     labs(x = "",
          y = "",
          fill = "log2(Fold_Change)")) +
    ggsave(paste(fig.path, "/DEGs/LFC-deseq_treatment_lfc-km_clustered_heatmap.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 20, 
           width = 4, limitsize = F, bg = "transparent")

idx <- which(df$gene_ids %in% row.names(logFC_P))
(cpm_hmap <- df[idx,] %>% 
group_by(gene_ids, genotype, treatment) %>%
summarise(counts = mean(counts)) %>%
ungroup(.) %>%
cbind.data.frame(., cluster_k = lfc_treatment$cluster_k[match(.$gene_ids, row.names(lfc_treatment))]) %>%
 mutate(gene_ids = factor(as.character(gene_ids), levels = gene_clusters_lfc),
   cluster_k = factor(cluster_k, levels = cluster_sorted)) %>% 
     ggplot(aes(x = genotype, 
                   y = gene_ids, 
                   fill = saturate(counts))) +
     geom_raster(alpha = 1) +
     facet_grid(cluster_k ~ treatment,
                switch = "both", space = "free", scale = "free",
                labeller = labeller(label.idx)) +
     scale_fill_gradient(na.value = "darkgrey", 
                         low = "black", 
                         high = "yellow") +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
           strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
           legend.text = element_text(size = 6)) +
     labs(x = "", y = "",
         fill = "log2[counts]")) +
    ggsave(paste(fig.path, "/EDA_heatmap_treatment_lfc_km-clustered.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 20, 
           width = 5, limitsize = F, bg = "transparent")


hmap_composite <- cowplot::plot_grid(
  strip_cluster + theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()), 
  cpm_hmap + theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()), 
  lfc_hmap + theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()), 
  nrow = 1, 
  ncol = 3, 
  axis = "tblr", 
  rel_widths = c(.75, 2, 1.5), 
  rel_height = c(1),
  greedy = FALSE
)

ggsave(hmap_composite, 
  file = paste(fig.path, "/DEGs/LFC_summary_km_clustered.png", sep = ""),
  width = 5, height = 25, 
  units = "in", device = "png", dpi = 600, limitsize = FALSE, bg = "transparent")


# Volcano Plot
lfc_cutoff <- c(-.5, .5)
FDR_cutoff <- c(-log10(5e-2), -log10(5e-5), -log10(5e-10), -log10(5e-15))

# ADD gene names for thos that are in cluster and in the outer
df_lfc %>%
filter(abs(logFC) > 0.05, genotype != "WT") %>%
mutate(gene_id = factor(gene_id, levels = gene_clusters_lfc),
    cluster_k = factor(cluster_k, levels = as.character(hc$order)),
    logFC = ifelse(abs(logFC)>10, 10*sign(logFC), logFC) , 
    mark = (abs(logFC) > 1 & PValue < 0.01),
    outer = (abs(logFC) > 1 & PValue < 1e-10)) %>%
    ggplot(aes(y = (-log10(PValue)), x = (logFC))) +
    geom_hline(yintercept = FDR_cutoff, lty = "solid", colour = "darkgrey", lwd = 0.8) +
    geom_vline(xintercept = lfc_cutoff, lty = "solid", colour = "darkgrey", lwd = 0.8) +
    geom_point(aes(alpha = mark, colour = outer, fill = mark, size = mark), shape = 23) +
    facet_grid(.~ genotype, scale = "fixed", space = "free", switch = "y") +
    scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "darkgrey"), guide = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = FALSE) +
    scale_size_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = FALSE) +
    scale_alpha_manual(values = c(`TRUE` = 0.8, `FALSE` = 0.4), guide = FALSE) +
    theme_AKB +
    theme(panel.spacing = unit(0.2, "lines"),
          strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
          axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
          axis.text.y = element_text(size = 16, angle = 30, vjust = 0, face = "bold")) +
    labs(x = "",
         y = "",
         fill = "log2(Fold_Change)") +
    ggsave(paste(fig.path, "/DEGs/LFC-deseq_treatment_lfc-volcanoplot.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 5, 
           width = 10, limitsize = F, bg = "transparent")


# Selected clusters only
sel_cluster <- c(10, 16, 5, 13)
clust_lfc_df <- df_lfc %>%
filter(cluster_k %in% sel_cluster) %>%
mutate(
  gene_id = factor(gene_id, levels = gene_clusters_lfc[(.$gene_id %in% gene_clusters_lfc)]),
  cluster_k = factor(cluster_k, levels = as.character(hc$order)),
  sig = ifelse(PValue < alpha & abs(logFC) > 1, "Significant", NA)) %>%
data.frame(.)

sig_genes <- as.character(unique(df_lfc$gene_id[which(df_lfc$PValue <= 0.05 & abs(df_lfc$logFC) >= 0.05)]))
lfc_mut <- logFC_P.annotated[sig_genes, !str_detect(colnames(logFC_P.annotated), "_WT_")]
ct <- cor.test(lfc_mut$treatment_myb_JA_logFC, lfc_mut$treatment_myc.tKO_JA_logFC)
d <- 1-cor(t(lfc_mut[,c("treatment_myb_JA_logFC", "treatment_myc.tKO_JA_logFC")]))
hc <- hclust(as.dist(d), "ward.D2")
sig_gene_clusters <- row.names(lfc_mut)[hc$order]

lfc_mut_clust <- lfc_mut %>% 
add_column(gene_ids = row.names(.), .before = 1) %>%
gather(key = "key", value= "lfc", convert = FALSE, treatment_myb_JA_logFC, treatment_myc.tKO_JA_logFC) %>%
separate(key, into = c("x", "Genotype", "x1", "x2"), convert = F, sep = "_") %>%
dplyr::select(-x, -x1, -x2) %>%
mutate(Genotype = factor(.$Genotype, levels = genotype$short[which(genotype$short %in% .$Genotype)])) %>%
data.frame(., stringsAsFactors = FALSE)


clust_lfc_df$ja <- str_detect(clust_lfc_df$gene_name, "JA|jasmonic|jasmonate")
clust_lfc_df$gsl <- str_detect(clust_lfc_df$gene_name, "gsl|glucosinolate")
clust_lfc_df$biosynthesis <- str_detect(clust_lfc_df$gene_name, "synth")
clust_lfc_df$trans <- str_detect(clust_lfc_df$gene_name, "transcript|transcribe|transcription")
clust_lfc_df$herbivore <- str_detect(clust_lfc_df$gene_name, "herbivore|bug|insect|feeding")
clust_lfc_df$pathogen <- str_detect(clust_lfc_df$gene_name, "pathogen|mutual|fungi")
clust_lfc_df$glucosidase <- str_detect(clust_lfc_df$gene_name, "glucosidase")
clust_lfc_df$ER <- str_detect(clust_lfc_df$gene_name, "ER|endoplasmic|reticulum")

(fig1a <- df %>% 
 group_by(gene_ids, genotype, treatment) %>%
 summarise(counts = mean(counts)) %>%
 ungroup(.) %>%
 cbind.data.frame(., cluster_k = lfc_treatment$cluster_k[match(.$gene_ids, row.names(lfc_treatment))]) %>%
 mutate(gene_ids = factor(as.character(gene_ids), levels = gene_clusters_lfc),
   cluster_k = factor(cluster_k, levels = as.character(cluster_sorted))) %>% 
 filter(cluster_k %in% sel_cluster) %>%
     ggplot(aes(x = genotype, 
                   y = gene_ids, 
                   fill = saturate(counts))) +
     geom_raster(alpha = 1) +
     facet_grid(cluster_k ~ treatment,
                switch = "both", space = "free", scale = "free",
                labeller = labeller(label.idx)) +
     scale_fill_gradient(low = "black", 
                         high = "yellow") +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
           strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
           legend.text = element_text(size = 6)) +
     labs(x = "", y = "",
         fill = "log2[counts]")) +
    ggsave(paste(fig.path, "/EDA_heatmap_treatment_lfc_sel_km-clustered.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 30, 
           width = 7, limitsize = F, bg = "transparent")

(fig1b <- clust_lfc_df %>%
  mutate(gene_ids = factor(as.character(gene_id), levels = gene_clusters_lfc),
   cluster_k = factor(cluster_k, levels = as.character(cluster_sorted))) %>% 
     ggplot(aes(x = genotype, y = gene_ids, 
                   fill = saturate(logFC))) +
     geom_raster(alpha = 1) +
     geom_tile(fill = NA, width = 0.95, height = 0.95, aes(colour = sig)) +
     facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
     scale_fill_gradient2(low = gradient.low, 
                          mid = gradient.mid,
                          high = gradient.high, 
                          na.value = gradient.na) +
     scale_colour_manual(values = c(`NS` = NA, `Significant` = "black")) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
           strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
           axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
           axis.text.y = element_blank()) +
     labs(x = "",
          y = "",
          fill = "log2(Fold_Change)")) +
    ggsave(paste(fig.path, "/DEGs/LFC-deseq_treatment_lfc-sel_km_clustered_heatmap.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 30, 
           width = 7, limitsize = F, bg = "transparent")

(fig2a <- df %>%
  filter(gene_ids %in% sig_genes) %>% 
 group_by(gene_ids, genotype, treatment) %>%
 summarise(counts = mean(counts)) %>%
 ungroup(.) %>%
 cbind.data.frame(., cluster_k = lfc_treatment$cluster_k[match(.$gene_ids, row.names(lfc_treatment))]) %>%
mutate(gene_ids = factor(as.character(gene_ids), levels = gene_clusters_lfc),
   cluster_k = factor(cluster_k, levels = as.character(cluster_sorted[cluster_sorted %in% .$cluster_k]))) %>% 
 # filter(cluster_k %in% sel_cluster) %>%
     ggplot(aes(x = genotype, 
                   y = gene_ids, 
                   fill = saturate(counts))) +
     geom_raster(alpha = 1) +
     facet_grid(cluster_k ~ treatment,
                switch = "both", space = "free", scale = "free",
                labeller = labeller(label.idx)) +
     scale_fill_gradient(low = "black", 
                         high = "yellow") +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
           strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
           legend.text = element_text(size = 6)) +
     labs(x = "", y = "",
         fill = "log2[counts]")) +
    ggsave(paste(fig.path, "/EDA_heatmap_treatment_lfc_sig_genes_km-clustered.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 7, 
           width = 3, limitsize = F, bg = "transparent")

(fig2b <- lfc_mut_clust %>%
    mutate(
      gene_ids = factor(gene_ids, levels = gene_clusters_lfc),
      cluster_k = factor(cluster_k, levels = as.character(cluster_sorted))) %>% 
     ggplot(aes(x = Genotype, y = gene_ids)) +
     geom_raster(alpha = 1, 
                   aes(fill = saturate(lfc))) +
     facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
     scale_fill_gradient2(low = gradient.low, 
                          mid = gradient.mid,
                          high = gradient.high, 
                          na.value = gradient.na) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
           strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
           axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
           axis.text.y = element_blank()) +
     labs(x = "",
          y = "",
          fill = "log2(Fold_Change)")) +
    ggsave(paste(fig.path, "/DEGs/LFC-deseq_treatment_lfc_sig_genes_km_clustered_heatmap.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 7, 
           width = 3, limitsize = F, bg = "transparent")

# go_term <-  unique(unlist(str_split(clust_lfc_df$go_id, ";")))
# term <- c("ja", "gsl", "biosynthesis", "trans", "herbivore", "pathogen", "glucosidase", "ER")
# plot_obj <- mclapply(term, function(x){

#       # idg <- str_detect(as.character(clust_lfc_df$go_id), pattern = x)
#       temp <- clust_lfc_df %>%
#       mutate(y = clust_lfc_df[,x]) %>%
#       filter(y == TRUE) %>%
#       data.frame(., stringsAsFactors = FALSE)

#       ida <- which(clust_lfc_df$genotype == "WT")
#       idb <- which(clust_lfc_df$genotype == "myc.tKO")
#       idc <- which(clust_lfc_df$genotype == "myb")

#       # mt_wt <- quantile(clust_lfc_df$logFC[ida],c(0.25,0.75))

#       # stat <- list(A_B = broom::tidy(cor.test(y = temp$logFC[idb], x = temp$logFC[ida], alternative = "two.sided")),
#       #             A_C =  broom::tidy(cor.test(y = temp$logFC[idc], x = temp$logFC[ida], alternative = "two.sided")),
#       #             B_C = broom::tidy(cor.test(y = temp$logFC[idb], x = temp$logFC[idc], alternative = "two.sided"))
#       #             )
      
#       # stat <- do.call(rbind.data.frame, stat) %>% 
#       # mutate(groups = names(stat), 
#       #   term = x) %>% 
#       # data.frame(.)

#       # message(print(stat))

#     (plt <- temp %>%
#         mutate(sig = (PValue < 0.01 & abs(logFC) > .5), 
#           aracyc = str_replace_all(aracyc, "^NA$", "Undetermined"),
#           symbol = str_replace_all(symbol, "^NA$", ""),
#           ) %>%
#       mutate(mark_text = ifelse(sig == TRUE, as.character(symbol), "")) %>%
#       arrange(desc(logFC)) %>%
#       ggplot(aes(x = as.factor(aracyc), 
#         y = logFC)
#       ) +
#       geom_hline(yintercept = c(-0.5, 0.5), lty = "solid", colour = "darkgrey", lwd = 0.8) +
#        geom_point(aes(alpha = sig, fill = sig, size = sig, colour = sig), shape = 23) +
#        ggrepel::geom_text_repel(aes(label = mark_text), size = 2) +
#        facet_grid(. ~ genotype, scale = "fixed", space = "free", switch = "x") +
#       scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "darkgrey"), guide = FALSE) +
#       scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = FALSE) +
#       scale_size_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = FALSE) +
#       scale_alpha_manual(values = c(`TRUE` = 0.8, `FALSE` = 0.4), guide = FALSE) +
#       # scale_shape_manual(values = c(`14` = 23, `15` = 24), guide = FALSE) +
#        theme_AKB +
#        theme(panel.spacing = unit(0.2, "lines"),
#              strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
#              axis.text.x = element_text(size = 6, angle = 0, vjust = 0, face = "bold"),
#              axis.text.y = element_text(size = 2, angle = 0, vjust = 0, face = "bold")) +
#        labs(x = "",
#             y = "") +
#        coord_flip()) +
#       ggsave(paste(fig.path, "/aracyc/term_", x,".png", sep = ""), 
#              dpi = 600, 
#              device = "png", 
#              units = "in",
#              height = 3, 
#              width = 6, limitsize = F, bg = "transparent")

#       obj <- list(plt = plt)
#       return(obj)
    

# }, mc.cores = 8)

# All aracyc terms
ida <- which(clust_lfc_df$genotype == "WT")
idb <- which(clust_lfc_df$genotype == "myc.tKO")
idc <- which(clust_lfc_df$genotype == "myb")

# mt_wt <- quantile(clust_lfc_df$logFC[ida],c(0.25,0.75))
clust_lfc_df %>%
      filter(aracyc != "NA", genotype != "WT") %>%
      mutate(sig = (PValue < 0.01 & abs(logFC) > .5),
        symbol = str_replace_all(symbol, "^NA$", ""),
        ) %>%
    # filter(aracyc != "Undetermined") %>%
    mutate(mark_text = ifelse(sig == TRUE, as.character(symbol), "")) %>%
    arrange(desc(logFC)) %>%
    ggplot(aes(x = as.factor(aracyc), 
      y = ifelse(abs(logFC) > 10, 10*sign(logFC), logFC))
    ) +
    geom_hline(yintercept = c(-0.5, 0.5), lty = "solid", colour = "darkgrey", lwd = 0.8) +
     geom_point(aes(alpha = sig, fill = sig, size = sig, colour = sig), shape = 23) +
     ggrepel::geom_text_repel(aes(label = mark_text), nudge_x = 0.05, size = 2) +
     facet_grid(. ~ genotype, scale = "fixed", space = "free", switch = "x") +
    scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "darkgrey"), guide = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = FALSE) +
    scale_size_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = FALSE) +
    scale_alpha_manual(values = c(`TRUE` = 0.8, `FALSE` = 0.4), guide = FALSE) +
    # scale_shape_manual(values = c(`14` = 23, `15` = 24), guide = FALSE) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
           strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
           axis.text.x = element_text(size = 12, angle = 0, vjust = 0, face = "bold"),
           axis.text.y = element_text(size = 1, angle = 0, vjust = 0, face = "bold")) +
     labs(x = "",
          y = "") +
     coord_flip() +
    ggsave(paste(fig.path, "/aracyc/all_term.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 10, 
           width = 15, limitsize = F, bg = "transparent")


lfc_mut_clust %>%
      filter(aracyc != "NA", Genotype != "WT") %>%
      mutate(sig = (abs(lfc) > 2),
        symbol = str_replace_all(symbol, "^NA$", ""),
        ) %>%
    # filter(aracyc != "Undetermined") %>%
    mutate(mark_text = ifelse(sig == TRUE, as.character(symbol), "")) %>%
    arrange(desc(lfc)) %>%
    ggplot(aes(x = as.factor(aracyc), 
      y = ifelse(abs(lfc) > 10, 10*sign(lfc), lfc))
    ) +
    geom_hline(yintercept = c(-0.5, 0.5), lty = "solid", colour = "darkgrey", lwd = 0.8) +
     geom_point(aes(alpha = sig, fill = sig, size = sig, colour = sig), shape = 23) +
     ggrepel::geom_text_repel(aes(label = mark_text), nudge_x = 0.05, size = 2) +
     facet_grid(. ~ Genotype, scale = "fixed", space = "free", switch = "x") +
    scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "darkgrey"), guide = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = FALSE) +
    scale_size_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = FALSE) +
    scale_alpha_manual(values = c(`TRUE` = 0.8, `FALSE` = 0.4), guide = FALSE) +
    # scale_shape_manual(values = c(`14` = 23, `15` = 24), guide = FALSE) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
           strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
           axis.text.x = element_text(size = 12, angle = 0, vjust = 0, face = "bold"),
           axis.text.y = element_text(size = 1, angle = 0, vjust = 0, face = "bold")) +
     labs(x = "",
          y = "") +
     coord_flip() +
    ggsave(paste(fig.path, "/aracyc/sig_genes_all_term.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 7, 
           width = 8, limitsize = F, bg = "transparent")

# Correlation between the significantly differing genes in mutants
lfc_mut %>%
filter(symbol != "NA") %>%
mutate(sig = ((abs(treatment_myb_JA_logFC) >= 1) & (abs(treatment_myc.tKO_JA_logFC) >= 1)),
  sigs = ((sign(treatment_myb_JA_logFC) == sign(treatment_myc.tKO_JA_logFC)) & 
    (abs(treatment_myb_JA_logFC) >= 2) & 
    (abs(treatment_myc.tKO_JA_logFC) >= 2))) %>%
 mutate(mark_text = ifelse(sigs == TRUE, replace_na(as.character(symbol), ""), "")) %>%
ggplot(aes(x = saturate(treatment_myb_JA_logFC), y = saturate(treatment_myc.tKO_JA_logFC))) +
geom_hline(yintercept = c(-1, 1), lty = "solid", colour = "darkgrey", lwd = 0.8) +
geom_vline(xintercept = c(-1, 1), lty = "solid", colour = "darkgrey", lwd = 0.8) +
geom_point(aes(alpha = sig, fill = sig, size = sig, colour = sig), shape = 23) +
ggrepel::geom_text_repel(aes(label = mark_text), nudge_x = 0.05, size = 2) +
scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "darkgrey"), guide = FALSE) +
scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = FALSE) +
scale_size_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = FALSE) +
scale_alpha_manual(values = c(`TRUE` = 0.8, `FALSE` = 0.4), guide = FALSE) +
theme_AKB +
theme(panel.spacing = unit(0.2, "lines"),
     strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
     axis.text.x = element_text(size = 10, angle = 0, vjust = 0, face = "bold"),
     axis.text.y = element_text(size = 10, angle = 0, vjust = 0, face = "bold")) +
labs(x = "",
    y = "") +
ggsave(paste(fig.path, "/DEGs/correlation_DEG_lfc_mutants.png", sep = ""), 
     dpi = 600, 
     device = "png", 
     units = "in",
     height = 4, 
     width = 4, limitsize = F, bg = "transparent")


# Strips for annotation
(strip_ja <- clust_lfc_df %>%
     ggplot(aes(x = "", y = gene_id, 
                   fill = ja)) +
     geom_raster(alpha = 1) +
     # geom_tile(fill = NA, width = 0.95, height = 0.95, aes(colour = sig)) +
     facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
     scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "red"), guide = FALSE) +
     # scale_colour_manual(values = c(`NA` = "white", `Significant` = "black"), guide = FALSE) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()) +
     labs(x = "",
          y = ""))

(strip_gsl <- clust_lfc_df %>%
     ggplot(aes(x = "", y = gene_id, 
                   fill = gsl)) +
     geom_raster(alpha = 1) +
     # geom_tile(fill = NA, width = 0.95, height = 0.95, aes(colour = sig)) +
     facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
     scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "red"), guide = FALSE) +
     # scale_colour_manual(values = c(`NA` = "white", `Significant` = "black"), guide = FALSE) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()) +
     labs(x = "",
          y = ""))


(strip_biosynthesis <- clust_lfc_df %>%
     ggplot(aes(x = "", y = gene_id, 
                   fill = biosynthesis)) +
     geom_raster(alpha = 1) +
     # geom_tile(fill = NA, width = 0.95, height = 0.95, aes(colour = sig)) +
     facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
     scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "red"), guide = FALSE) +
     # scale_colour_manual(values = c(`NA` = "white", `Significant` = "black"), guide = FALSE) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()) +
     labs(x = "",
          y = ""))


(strip_trans <- clust_lfc_df %>%
     ggplot(aes(x = "", y = gene_id, 
                   fill = trans)) +
     geom_raster(alpha = 1) +
     # geom_tile(fill = NA, width = 0.95, height = 0.95, aes(colour = sig)) +
     facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
     scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "red"), guide = FALSE) +
     # scale_colour_manual(values = c(`NA` = "white", `Significant` = "black"), guide = FALSE) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()) +
     labs(x = "",
          y = ""))

(strip_herbivore <- clust_lfc_df %>%
     ggplot(aes(x = "", y = gene_id, 
                   fill = herbivore)) +
     geom_raster(alpha = 1) +
     # geom_tile(fill = NA, width = 0.95, height = 0.95, aes(colour = sig)) +
     facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
     scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "red"), guide = FALSE) +
     # scale_colour_manual(values = c(`NA` = "white", `Significant` = "black"), guide = FALSE) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()) +
     labs(x = "",
          y = ""))

(strip_pathogen <- clust_lfc_df %>%
     ggplot(aes(x = "", y = gene_id, 
                   fill = pathogen)) +
     geom_raster(alpha = 1) +
     # geom_tile(fill = NA, width = 0.95, height = 0.95, aes(colour = sig)) +
     facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
     scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "red"), guide = FALSE) +
     # scale_colour_manual(values = c(`NA` = "white", `Significant` = "black"), guide = FALSE) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()) +
     labs(x = "",
          y = ""))

(strip_glucosidase <- clust_lfc_df %>%
     ggplot(aes(x = "", y = gene_id, 
                   fill = glucosidase)) +
     geom_raster(alpha = 1) +
     # geom_tile(fill = NA, width = 0.95, height = 0.95, aes(colour = sig)) +
     facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
     scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "red"), guide = FALSE) +
     # scale_colour_manual(values = c(`NA` = "white", `Significant` = "black"), guide = FALSE) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()) +
     labs(x = "",
          y = ""))

(strip_ER <- clust_lfc_df %>%
     ggplot(aes(x = "", y = gene_id, 
                   fill = ER)) +
     geom_raster(alpha = 1) +
     # geom_tile(fill = NA, width = 0.95, height = 0.95, aes(colour = sig)) +
     facet_grid(cluster_k ~., scale = "free", space = "free", switch = "y") +
     scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "red"), guide = FALSE) +
     # scale_colour_manual(values = c(`NA` = "white", `Significant` = "black"), guide = FALSE) +
     theme_AKB +
     theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()) +
     labs(x = "",
          y = ""))

hmap_composite <- cowplot::plot_grid(
  fig1a + theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()), 
  fig1b + theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()), 
  strip_ja,
  strip_gsl,
  strip_biosynthesis,
  strip_trans,
  strip_herbivore,
  strip_pathogen,
  strip_glucosidase,
  strip_ER,
  nrow = 1, 
  ncol = 10, 
  axis = "tblr", 
  rel_widths = c(0.65, .5, .3, .3, .3, .3, .3, .3, .3, .3), 
  rel_height = c(1),
  greedy = FALSE
)


ggsave(hmap_composite, 
  file = paste(fig.path, "/DEGs/LFC_summary_annotated_sel_km_clustered.png", sep = ""),
  width = 7, height = 28, 
  units = "in", device = "png", dpi = 600, limitsize = FALSE, bg = "transparent")

# Composite figures of the significant genes only
hmap_composite <- cowplot::plot_grid(
  fig2a + theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()), 
  fig2b + theme(panel.spacing = unit(0.2, "lines"),
    legend.position = "none", 
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()), 
  nrow = 1, 
  ncol = 2, 
  axis = "tblr", 
  rel_widths = c(0.65, .5), 
  rel_height = c(1),
  greedy = FALSE
)

ggsave(hmap_composite, 
  file = paste(fig.path, "/DEGs/LFC_summary_annotated_sig_genes_km_clustered.png", sep = ""),
  width = 2, height = 6, 
  units = "in", device = "png", dpi = 600, limitsize = FALSE, bg = "transparent")


# Differential expression within the cluster
lfc_cutoff <- c(-0.5, 0.5)
FDR_cutoff <- c(-log10(5e-2), -log10(5e-5))
# wt_mark <- range(clust_lfc_df$logFC[clust_lfc_df$genotype == "WT" & clust_lfc_df$PValue < 0.05 & abs(clust_lfc_df$logFC) > 0.5])

# ADD gene names for thos that are in cluster and in the outer
clust_lfc_df %>%
# filter(genotype != "WT") %>%
mutate(
    symbol = str_replace_all(symbol, "^NA$", ""),
    gene_id = factor(gene_id, levels = gene_clusters_lfc),
    cluster_k = factor(cluster_k, levels = as.character(hc$order)),
    logFC = ifelse(abs(logFC)>10, 10*sign(logFC), logFC) , 
    mark = (abs(logFC) > 0.5 & PValue < 0.01),
    outer = (gene_id %in% target)) %>%
mutate(mark_text = ifelse(outer == TRUE & !is.na(symbol), as.character(symbol), "")) %>%
    ggplot(aes(y = (-log10(PValue)), x = (logFC))) +
    geom_hline(yintercept = FDR_cutoff, lty = "solid", colour = "darkgrey", lwd = 0.8) +
    geom_vline(xintercept = c(lfc_cutoff), lty = "solid", colour = "darkgrey", lwd = 0.8) +
    geom_point(aes(alpha = mark, colour = outer, fill = mark, size = mark), shape = 23) +
    ggrepel::geom_text_repel(aes(label = mark_text), size = 2) +
    facet_grid(.~ genotype, scale = "fixed", space = "free", switch = "y") +
    scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "darkgrey"), guide = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = FALSE) +
    scale_size_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = FALSE) +
    scale_alpha_manual(values = c(`TRUE` = 0.8, `FALSE` = 0.4), guide = FALSE) +
    theme_AKB +
    theme(panel.spacing = unit(0.2, "lines"),
          strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
          axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
          axis.text.y = element_text(size = 16, angle = 30, vjust = 0, face = "bold")) +
    labs(x = "",
         y = "",
         fill = "log2(Fold_Change)") +
    ggsave(paste(fig.path, "/DEGs/LFC-deseq_treatment_lfc_sel_km-clustered-volcanoplot.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 5, 
           width = 10, limitsize = F, bg = "transparent")

# Plot individual clusters of gene



# Iterate over the cluster14 and cluster15 find the names of the genes and ARACYC function
# Cumulative gene expression of these genes with respect to ARACYC terms 
# Mark GENE names if they are significant

message("Spilling Venn Diagrams for differentially regulated genes ...")
# Plot Venn diagram for differential genes
idx <- str_detect(colnames(DEGs), "^treatment")
df_DEG <- DEGs[,idx]

regulated <- list(
    wt_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "WT_JA")] != 0)], 
    myb_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myb_JA")] != 0)], 
    myc_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myc.tKO_JA")] != 0)] 
)

# DEGs
VennDiagram::venn.diagram(x = regulated,
                          category.names = label.idx,
                          filename = paste(fig.path, "/DEGs/DEGs-all_deseq_treatment_VenD.png", sep = ""),
                          imagetype = "png",
                          sigdigs = 2,
                          hyper.test = TRUE, 
                          lower.tail = TRUE,
                          output = T,
                          height = 1000,
                          width = 1000,
                          resolution = 600,
                          lwd = 3,
                          lty = "blank",
                          fill = genotype$col.idx,
                          cex = 0.8,
                          fontface = "bold",
                          fontfamily = "sans",
                          cat.cex = 0.4,
                          cat.fontface = "bold",
                          cat.default.pos = "outer", 
                          main = paste0("Differentially Expressed Genes(logFC threshold = ", 
                                        FC_threshold,"; FDR <", alpha,") - Effect of treatment"), 
                          #print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)
upregulated <- list(
    wt_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "WT_JA")] == 1)], 
    myb_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myb_JA")] == 1)], 
    myc_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myc.tKO_JA")] == 1)]
    )

# Upregulated
VennDiagram::venn.diagram(x = upregulated,
                          category.names = label.idx,
                          filename = paste(fig.path, "/DEGs/DEGs-up_deseq_treatment_VenD.png", sep = ""),
                          imagetype = "png",
                          sigdigs = 2,
                          hyper.test = TRUE, 
                          lower.tail = TRUE,
                          output = T,
                          height = 1000,
                          width = 1000,
                          resolution = 600,
                          lwd = 3,
                          lty = "blank",
                          fill = genotype$col.idx,
                          cex = 0.8,
                          fontface = "bold",
                          fontfamily = "sans",
                          cat.cex = 0.4,
                          cat.fontface = "bold",
                          cat.default.pos = "outer", 
                          main = paste0("Upregulated Genes(logFC threshold = ", 
                                        FC_threshold,"; FDR <", alpha,") - Effect of treatment"), 
                          #print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)
downregulated <- list(
    wt_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "WT_JA")] == -1)], 
    myb_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myb_JA")] == -1)], 
    myc_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myc.tKO_JA")] == -1)] 
)

# Downregulated
VennDiagram::venn.diagram(x = downregulated,
                          category.names = label.idx,
                          filename = paste(fig.path, "/DEGs/DEGs-down_deseq_treatment_VenD.png", sep = ""),
                          imagetype = "png",
                          sigdigs = 2,
                          hyper.test = TRUE, 
                          lower.tail = TRUE,
                          output = T,
                          height = 1000,
                          width = 1000,
                          resolution = 600,
                          lwd = 3,
                          lty = "blank",
                          fill = genotype$col.idx,
                          cex = 0.8,
                          fontface = "bold",
                          fontfamily = "sans",
                          cat.cex = 0.4,
                          cat.fontface = "bold",
                          cat.default.pos = "outer", 
                          main = paste0("Downregulated Genes(logFC threshold = ", 
                                        FC_threshold,"; FDR <", alpha,") - Effect of treatment"), 
                          #print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)
message("DONE")
# save list of genes in an object for GSET analysis
save(list = c("regulated", "downregulated", "upregulated"),
     file = paste(out.path, "/DEG-list_deseq-treatment.Rdata", sep = ""))

# END OF SCRIPT
sessionInfo()
