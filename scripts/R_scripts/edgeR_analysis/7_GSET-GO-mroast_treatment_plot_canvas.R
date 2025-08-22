rm(list = ls())

# NOTE::
# Plot results from LFC GSET ANALYSIS
# Analysis from MROAST
# Change lfc terms to FDR

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
load(file = paste(out.path, "GSET-GO_edger_objects-treatment.plotting.Rdata", sep = ""))

# Make wide data
idx <- str_detect(mroast.GO.long$contrast, "^treatment") 
idy <- str_detect(colnames(mroast.GO.long), "^FDR|go_id|contrast|aracyc|gene_name")
go_p <- mroast.GO.long[idx,idy]
go_p <- go_p %>% 
            spread(key = "contrast", value = "FDR", convert = FALSE) %>%
            data.frame(.)
row.names(go_p) <- go_p$go_id
idx <- str_detect(colnames(go_p), "^treatment") 
d <- as.dist(1 - cor(t(go_p[,idx])))
clust_go_fdr <- hclust(as.dist(d), method = "ward.D")
go_clusters_fdr <- row.names(go_p[,idx])[clust_go_fdr$order]

idx <- str_detect(colnames(go_p), "WT_JA") 
idy <- order(go_p[,idx])
gene_sorted_wt <- row.names(go_p)[idy]

df_path <- go_p %>% 
            dplyr::select(go_id, aracyc) %>%
            data.frame(., row.names = .$go_id)

path.sorted <- sort(names(path.idx))
# Index out the pathway definitions
for(i in 1:length(names(path.idx))){
    x <- path.idx[[i]]
    nam <- names(path.idx)[i]
    char <- ifelse(length(x) > 1, paste(x, collapse = "|"), x)
    df_path[,nam] = 0
    idx <- str_detect(df_path$aracyc, pattern = char)
    df_path[idx, nam] <- 1
}

# Differential Gene Data carpenting
idx <- str_detect(mroast.GO.long$contrast, "^treatment")
go_p <- mroast.GO.long[idx,]
df_fdr <- go_p %>% 
            separate(contrast, into = c("x1", "genotype", "x2"), 
                     convert = FALSE, remove = FALSE, sep = "_") %>%
            dplyr::select(-x1,-x2) %>%
            mutate(go_id = factor(go_id, levels = go_clusters_fdr)) %>%
            data.frame(.)


# Make factors
df_fdr$genotype <- factor(df_fdr$genotype, levels= genotype$short)
#df_fdr$treatment <- factor(df_fdr$treatment, levels= treatment$name)
df_fdr$Direction <- as.numeric(ifelse(df_fdr$Direction == "Up", 1, -1))
df_fdr$sig <- as.factor(ifelse(df_fdr$FDR < alpha, 1 * df_fdr$Direction, 0))
df_fdr[,path.sorted] <- df_path[match(df_fdr$go_id, row.names(df_path)), path.sorted]
df_fdr <- df_fdr %>% gather(key = "definition", 
                            val = "vals", 
                            convert = FALSE, 
                            path.sorted) %>%
            dplyr::filter(vals != 0) %>%
            dplyr::select(-vals) %>%
            mutate(definition = factor(definition, levels = path.sorted)) %>%
            data.frame(.)

message("Spilling GO heatmap gold standards ...")
# Those which are significant
df_fdr %>% ggplot(aes(x = genotype, y = go_id, fill = sig)) +
    geom_tile(alpha = 1, width = 0.6) +
    facet_grid(definition ~. , switch = "y", scales = "free", space = "free") +
    scale_fill_manual(values = c(`0` = gradient.na,
                                 `1` = gradient.sig.high, 
                                 `-1` = gradient.sig.low),
                     labels = c(`0` = "Not Significant",
                                 `1` = "Upregulated", 
                                `-1` = "Downregulated")) +
    theme_AKB +
    theme(axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
          axis.text.y = element_blank(),
          legend.text = element_text(size = 8), 
          strip.text.y = element_text(size = 8, angle = 180, hjust = 0.5, vjust = 0.5)) +
    labs(x = "",
         y = "",
         fill = paste("vs DW Treatment (FDR",alpha,")", sep = "")) +
    ggsave(paste(fig.path, 
                 "/GSET/GSET-GO-mroast-FDR[", alpha,"]-edger_treatment_GO-clustered_heatmap.png", 
                 sep = ""), 
           dpi = 300, 
           device = "png", 
           height = 30, 
           width = 10, limitsize = F)

df_fdr %>% mutate(go_id = factor(go_id, levels = gene_sorted_wt)) %>%
    ggplot(aes(x = genotype, y = go_id, fill = sig)) +
    geom_tile(alpha = 1, width = 0.6) +
    facet_grid(definition ~. , switch = "y", scales = "free", space = "free") +
    scale_fill_manual(values = c(`0` = gradient.na,
                                 `1` = gradient.sig.high, 
                                 `-1` = gradient.sig.low),
                     labels = c(`0` = "Not Significant",
                                 `1` = "Upregulated", 
                                `-1` = "Downregulated")) +
    theme_AKB +
    theme(axis.text.x = element_text(size = 16, angle = 30, vjust = 0, face = "bold"),
          axis.text.y = element_blank(),
          legend.text = element_text(size = 6),
          strip.text.y = element_text(size = 8, angle = 180, hjust = 0.5, vjust = 0.5)) +
    labs(x = "",
         y = "",
         fill =paste("vs DW Treatment (FDR",alpha,")", sep = "")) +
    ggsave(paste(fig.path, 
                 "/GSET/GSET-GO-mroast-FDR[", alpha,"]-edger_treatment_GO-WT_sorted_heatmap.png", 
                 sep = ""), 
           dpi = 300, 
           device = "png", 
           height = 30, 
           width = 10, limitsize = F)
    
message("Spilling Venn Diagrams for differentially regulated genes ...")
# Plot Venn diagram for differential genes
# Make a wide matrix of the df
df_DGO <- go_p %>% 
            dplyr::select(go_id, contrast, Direction, FDR) %>% 
            mutate(sig = as.numeric(ifelse(Direction == "Up", 1 * FDR, -1 * FDR))) %>%
            mutate(sig = ifelse(abs(sig) < alpha, sign(sig), 0)) %>%
            dplyr::select(go_id, contrast, sig) %>%
            spread(key = "contrast", value = "sig", convert = FALSE) %>%
            data.frame(.,row.names = .$go_id)

regulated <- list(
    wt_ja = row.names(df_DGO)[which(df_DGO[,str_detect(colnames(df_DGO), "WT_JA")] != 0)], 
    myb_ja = row.names(df_DGO)[which(df_DGO[,str_detect(colnames(df_DGO), "myb_JA")] != 0)], 
    myc_ja = row.names(df_DGO)[which(df_DGO[,str_detect(colnames(df_DGO), "myc.tKO_JA")] != 0)] 
)

# DEGs
VennDiagram::venn.diagram(x = regulated,
                          category.names = label.idx,
                          filename = paste(fig.path, 
                                           "/GSET/GSET-GO-mroast_all_treatment_VenD.png", 
                                           sep = ""),
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
                          main = paste("Differentially regulated GO:TERMS - 
Effect of treatment (FDR < ", alpha, ")", sep = ""), 
                          print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)

upregulated <- list(
    wt_ja = row.names(df_DGO)[which(df_DGO[,str_detect(colnames(df_DGO), "WT_JA")] == 1)], 
    myb_ja = row.names(df_DGO)[which(df_DGO[,str_detect(colnames(df_DGO), "myb_JA")] == 1)], 
    myc_ja = row.names(df_DGO)[which(df_DGO[,str_detect(colnames(df_DGO), "myc.tKO_JA")] == 1)] 
)

# DEGs
VennDiagram::venn.diagram(x = upregulated,
                          category.names = label.idx,
                          filename = paste(fig.path, 
                                           "/GSET/GSET-GO-mroast_up_treatment_VenD.png", 
                                           sep = ""),
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
                          main = paste("Differentially upregulated GO:TERMS - 
Effect of treatment (FDR < ", alpha, ")", sep = ""), 
                          print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)

downregulated <- list(
    wt_ja = row.names(df_DGO)[which(df_DGO[,str_detect(colnames(df_DGO), "WT_JA")] == -1)], 
    myb_ja = row.names(df_DGO)[which(df_DGO[,str_detect(colnames(df_DGO), "myb_JA")] == -1)], 
    myc_ja = row.names(df_DGO)[which(df_DGO[,str_detect(colnames(df_DGO), "myc.tKO_JA")] == -1)] 
)

# DEGs
VennDiagram::venn.diagram(x = downregulated,
                          category.names = label.idx,
                          filename = paste(fig.path, 
                                           "/GSET/GSET-GO-mroast_down_treatment_VenD.png", 
                                           sep = ""),
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
                          main = paste("Differentially downregulated GO:TERMS - 
Effect of treatment (FDR < ", alpha, ")", sep = ""), 
                          print.mode = "percent",
                          main.cex = 0.4
                          # cat.pos = c(-27, 27, 135),
                          # cat.dist = c(0.055, 0.055, 0.085),
                          # cat.fontfamily = "sans",
                          # rotation = 1
)

# save list of genes in an object for GSET analysis
save(list = c("regulated", "downregulated", "upregulated"),
     file = paste(out.path, "/GO-list-mroast_edger-treatment.Rdata", sep = ""))

# END OF SCRIPT
sessionInfo()
