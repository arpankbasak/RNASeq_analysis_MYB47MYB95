rm(list = ls())

options(warn = 1,
       mc.cores = 8)

pkgs <- c("tidyverse", "VennDiagram")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
load(paste(out.path, "/edger_objects.plotting.Rdata", sep = ""))

# Cluster genes on the basis of LFC- expression
idx <- str_detect(colnames(logFC_P), "_logFC$")
d <- 1 - cor(t(logFC_P[,idx]))
clust_gene_lfc <- hclust(as.dist(d), method = "ward.D")
gene_clusters_lfc <- row.names(logFC_P[,idx])[clust_gene_lfc$order]
# save(list = c("clust_gene_lfc", "gene_clusters_lfc"),
#      file = paste(out.path, "/LFC-cluster_edger-treatment.Rdata", sep = ""))

# Differential Gene Data carpenting
idx <- str_detect(colnames(logFC_P), "^genotype")
df_lfc <- logFC_P[,idx] %>% 
            add_column(gene_id = row.names(logFC_P[,idx]), .before = 1) %>%
            gather(key = "key", value = "val", convert = FALSE, -gene_id) %>%
            separate(key, into = c("comparison", "x1", "genotype", "treatment", "col"), convert = FALSE, sep = "_") %>%
            spread(key = "col", value = "val", convert = FALSE) %>%
            dplyr::select(-x1) %>%
            mutate(gene_id = factor(gene_id, levels = gene_clusters_lfc)) %>%
            data.frame(.)

# Make factors
df_lfc$genotype <- factor(df_lfc$genotype, levels= genotype$short)
df_lfc$sig <- as.factor(ifelse(df_lfc$PValue < alpha, 1 * sign(df_lfc$logFC), 0))


# Volcano plot
lfc_cutoff <- c(-.5, .5)
FDR_cutoff <- c(-log10(5e-2), -log10(5e-5), -log10(5e-10), -log10(5e-15))

df_lfc %>%
filter(abs(logFC) > 0.05) %>%
mutate(
    gene_id = factor(gene_id, levels = gene_clusters_lfc),
    # cluster_k = factor(cluster_k, levels = as.character(hc$order)),
    logFC = ifelse(abs(logFC)>10, 10*sign(logFC), logFC) , 
    mark = (abs(logFC) > 1 & PValue < 0.01),
    outer = (abs(logFC) > 1 & PValue < 1e-10),
    comparison = as.factor(comparison)) %>%
    ggplot(aes(y = (-log10(PValue)), x = (logFC))) +
    geom_hline(yintercept = FDR_cutoff, lty = "solid", colour = "darkgrey", lwd = 0.8) +
    geom_vline(xintercept = lfc_cutoff, lty = "solid", colour = "darkgrey", lwd = 0.8) +
    geom_point(aes(alpha = mark, colour = outer, fill = mark, size = mark), shape = 23) +
    facet_grid(.~ comparison + genotype, scale = "fixed", space = "free", switch = "y") +
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
    ggsave(paste(fig.path, "/DEGs/LFC-edger_lfc-volcanoplot.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           units = "in",
           height = 5, 
           width = 10, limitsize = F, bg = "transparent")

#

message("Spilling LFC heatmap ...")
# PLot Canvas for LFC
df_lfc %>% ggplot(aes(x = genotype, y = gene_id, fill = saturate(logFC))) +
    geom_tile(alpha = 1, width = 0.6) +
    scale_fill_gradient2(low = gradient.low, 
                         mid = gradient.mid,
                         high = gradient.high, 
                         na.value = gradient.na) +
    facet_grid(.~ treatment, scales = "free", space = "free") + 
    theme_AKB +
    theme(axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
          axis.text.y = element_blank(),
          legend.text = element_text(size = 6)) +
    labs(x = "",
         y = "",
         fill = "log2(Fold_Change)") +
    ggsave(paste(fig.path, "/DEGs/LFC-edger_overall_lfc-clustered_heatmap.png", sep = ""), 
           dpi = 300, 
           device = "png", 
           height = 30, 
           width = 10, limitsize = F)

message("Spilling LFC heatmap gold standards ...")
# Those which are significant
df_lfc %>% ggplot(aes(x = genotype, y = gene_id, fill = sig)) +
    geom_tile(alpha = 1, width = 0.6) +
    facet_grid(.~ treatment, scales = "free", space = "free") + 
    scale_fill_manual(values = c(`0` = gradient.na,
                                 `1` = gradient.sig.high, 
                                 `-1` = gradient.sig.low),
                     labels = c(`0` = "Not Significant",
                                 `1` = "Upregulated", 
                                 `-1` = "Downregulated")) +
    theme_AKB +
    theme(axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
          axis.text.y = element_blank(),
          legend.text = element_text(size = 6)) +
    labs(x = "",
         y = "",
         fill = paste("FDR < ", alpha, sep = "")) +
    ggsave(paste(fig.path, "/DEGs/LFC-sig[", alpha,"]-edger_overall_lfc-clustered_heatmap.png", sep = ""), 
           dpi = 300, 
           device = "png", 
           height = 30, 
           width = 10, limitsize = F)

df_lfc$high_sig <- ifelse(-log10(df_lfc$PValue) > 15, df_lfc$gene_id, NA)

message("Spilling Venn Diagrams for differentially regulated genes ...")
# Plot Venn diagram for differential genes
idx <- str_detect(colnames(DEGs), "^overall")
df_DEG <- DEGs[,idx]

regulated_ja <- list(
    wt_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "WT_JA")] != 0)], 
    myb_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myb_JA")] != 0)], 
    myc_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myc.tKO_JA")] != 0)] 
)
regulated_dw <- list(
    wt_dw = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "WT_DW")] != 0)], 
    myb_dw = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myb_DW")] != 0)], 
    myc_dw = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myc.tKO_DW")] != 0)] 
)

# DEGs
VennDiagram::venn.diagram(x = regulated_ja,
                          category.names = label.idx,
                          filename = paste(fig.path, "/DEGs/DEGs-all_edgeR_overall_JA_VenD.png", sep = ""),
                          imagetype = "png",
                          sigdigs = 2,
                          hyper.test = TRUE, 
                          lower.tail = TRUE,
                          output = T,
                          height = 1000,
                          width = 1000,
                          resolution = 300,
                          lwd = 3,
                          lty = "blank",
                          fill = genotype$col.idx,
                          cex = 0.2,
                          fontface = "bold",
                          fontfamily = "sans",
                          cat.cex = 0.4,
                          cat.fontface = "bold",
                          cat.default.pos = "outer", 
                          main = "Differentially Expressed Genes - Effect of JA", 
                          #print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)
VennDiagram::venn.diagram(x = regulated_dw,
                          category.names = label.idx,
                          filename = paste(fig.path, "/DEGs/DEGs-all_edgeR_overall_DW_VenD.png", sep = ""),
                          imagetype = "png",
                          sigdigs = 2,
                          hyper.test = TRUE, 
                          lower.tail = TRUE,
                          output = T,
                          height = 1000,
                          width = 1000,
                          resolution = 300,
                          lwd = 3,
                          lty = "blank",
                          fill = genotype$col.idx,
                          cex = 0.2,
                          fontface = "bold",
                          fontfamily = "sans",
                          cat.cex = 0.4,
                          cat.fontface = "bold",
                          cat.default.pos = "outer", 
                          main = "Differentially Expressed Genes - Effect of DW", 
                          #print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)
upregulated_ja <- list(
    wt_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "WT_JA")] == 1)], 
    myb_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myb_JA")] == 1)], 
    myc_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myc.tKO_JA")] == 1)]
    )
upregulated_dw <- list(
    wt_dw = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "WT_DW")] == 1)], 
    myb_dw = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myb_DW")] == 1)], 
    myc_dw = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myc.tKO_DW")] == 1)]
    )
# Upregulated
VennDiagram::venn.diagram(x = upregulated_ja,
                          category.names = label.idx,
                          filename = paste(fig.path, "/DEGs/DEGs-up_edgeR_overall_JA_VenD.png", sep = ""),
                          imagetype = "png",
                          sigdigs = 2,
                          hyper.test = TRUE, 
                          lower.tail = TRUE,
                          output = T,
                          height = 1000,
                          width = 1000,
                          resolution = 300,
                          lwd = 3,
                          lty = "blank",
                          fill = genotype$col.idx,
                          cex = 0.2,
                          fontface = "bold",
                          fontfamily = "sans",
                          cat.cex = 0.4,
                          cat.fontface = "bold",
                          cat.default.pos = "outer", 
                          main = "Upregulated Genes - Effect of JA", 
                          #print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)
VennDiagram::venn.diagram(x = upregulated_dw,
                          category.names = label.idx,
                          filename = paste(fig.path, "/DEGs/DEGs-up_edgeR_overall_JA_VenD.png", sep = ""),
                          imagetype = "png",
                          sigdigs = 2,
                          hyper.test = TRUE, 
                          lower.tail = TRUE,
                          output = T,
                          height = 1000,
                          width = 1000,
                          resolution = 300,
                          lwd = 3,
                          lty = "blank",
                          fill = genotype$col.idx,
                          cex = 0.2,
                          fontface = "bold",
                          fontfamily = "sans",
                          cat.cex = 0.4,
                          cat.fontface = "bold",
                          cat.default.pos = "outer", 
                          main = "Upregulated Genes - Effect of DW", 
                          #print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)
downregulated_ja <- list(
    wt_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "WT_JA")] == -1)], 
    myb_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myb_JA")] == -1)], 
    myc_ja = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myc.tKO_JA")] == -1)] 
)
downregulated_dw <- list(
    wt_dw = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "WT_DW")] == -1)], 
    myb_dw = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myb_DW")] == -1)], 
    myc_dw = row.names(df_DEG)[which(df_DEG[,str_detect(colnames(df_DEG), "myc.tKO_DW")] == -1)] 
)

# Downregulated
VennDiagram::venn.diagram(x = downregulated_ja,
                          category.names = label.idx,
                          filename = paste(fig.path, "/DEGs/DEGs-down_edgeR_overall_JA_VenD.png", sep = ""),
                          imagetype = "png",
                          sigdigs = 2,
                          hyper.test = TRUE, 
                          lower.tail = TRUE,
                          output = T,
                          height = 1000,
                          width = 1000,
                          resolution = 300,
                          lwd = 3,
                          lty = "blank",
                          fill = genotype$col.idx,
                          cex = 0.2,
                          fontface = "bold",
                          fontfamily = "sans",
                          cat.cex = 0.4,
                          cat.fontface = "bold",
                          cat.default.pos = "outer", 
                          main = "Downregulated Genes - Effect of JA", 
                          #print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)
VennDiagram::venn.diagram(x = downregulated_dw,
                          category.names = label.idx,
                          filename = paste(fig.path, "/DEGs/DEGs-down_edgeR_overall_JA_VenD.png", sep = ""),
                          imagetype = "png",
                          sigdigs = 2,
                          hyper.test = TRUE, 
                          lower.tail = TRUE,
                          output = T,
                          height = 1000,
                          width = 1000,
                          resolution = 300,
                          lwd = 3,
                          lty = "blank",
                          fill = genotype$col.idx,
                          cex = 0.2,
                          fontface = "bold",
                          fontfamily = "sans",
                          cat.cex = 0.4,
                          cat.fontface = "bold",
                          cat.default.pos = "outer", 
                          main = "Downregulated Genes - Effect of DW", 
                          #print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)
message("DONE")
# save list of genes in an object for GSET analysis
save(list = c("regulated_dw", "regulated_ja", "downregulated_ja",  "downregulated_dw", "upregulated_ja", "upregulated_dw"),
     file = paste(out.path, "/DEG-list_edger-overall.Rdata", sep = ""))

# END OF SCRIPT
sessionInfo()
