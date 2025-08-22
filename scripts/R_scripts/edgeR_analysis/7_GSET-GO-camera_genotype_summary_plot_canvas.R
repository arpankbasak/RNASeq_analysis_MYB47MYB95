# Script for individual analysis of genes corresponding to pathway
# @Arpan Kumar Basak
# TODO::Make manhattan plot for each contrast and specify the significan genes at threshold

# Cleanup
rm(list = ls())
# Manhattan plot

pkgs <- c("tidyverse", "limma", "parallel")
lapply(pkgs, require, character.only = T)
options(set.cores = detectCores()/2)

path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

options(set.cores = detectCores()/2)

# Loading dependencies
load(paste(out.path, "edger_objects.plotting.Rdata", sep = "")) #LFC data
edata.cpm <- read.table(file = paste(stat.path, "/log2cpm.edger.txt", sep = ""), sep = "\t", header = TRUE) # ExpressionData # AnnotationData
load(file = paste(out.path, "GSET-GO_edger_objects-treatment.plotting.Rdata", sep = "")) # GO-Data-camera

# Cluster genes on the basis of LFC- expression -- use logFC_P.wide
logFC_P <- logFC_P
idx <- str_detect(colnames(logFC_P), "^genotype") & str_detect(colnames(logFC_P), "_logFC$")
d <- dist(logFC_P[,idx])
clust_gene_lfc <- hclust(d, method = "ward.D")
gene_clusters_lfc <- row.names(logFC_P[,idx])[clust_gene_lfc$order]

# Cluster ontologies on the basis of GO -- GO
idx <- str_detect(camera.GO.long$contrast, "^genotype") 
idy <- str_detect(colnames(camera.GO.long), "^FDR|go_id|contrast|aracyc|gene_name")
go_p <- camera.GO.long[idx,idy]
go_p <- go_p %>% 
            spread(key = "contrast", value = "FDR", convert = FALSE) %>%
            data.frame(.)
row.names(go_p) <- go_p$go_id
idx <- str_detect(colnames(go_p), "^genotype") 
d <- as.dist(1 - cor(t(go_p[,idx])))
clust_go_fdr <- hclust(as.dist(d), method = "ward.D")
go_clusters_fdr <- row.names(go_p[,idx])[clust_go_fdr$order]

# Significant DEGs
idx <- str_detect(logFC_P.long$contrast, "^genotype")
logFC_P <- logFC_P.long[idx,] %>% 
            mutate(gene_ids = factor(gene_ids, levels = gene_clusters_lfc)) %>%
            separate(contrast, into = c("x1", "x2", "genotype", "treatment"),
                                          sep = "_", convert = FALSE, remove = FALSE) %>%
            dplyr::select(-x1,-x2) %>% 
            data.frame(., stringsAsFactors = FALSE)

gene_sig <- as.character(unique(logFC_P$gene_ids[logFC_P$PValue < alpha]))

# Significant GO terms
idx <- str_detect(camera.GO.long$contrast, "^genotype")
go_p <- camera.GO.long[idx,] %>% 
            mutate(go_id = factor(go_id, levels = go_clusters_fdr)) %>%
            separate(contrast, into = c("x1", "x2", "genotype", "treatment"),
                                          sep = "_", convert = FALSE, remove = FALSE) %>%
            mutate(idx = paste(go_id, contrast, sep = "_")) %>%
            dplyr::select(-x1,-x2) %>% data.frame(., stringsAsFactors = FALSE)

go_sig <- as.character(unique(go_p$go_id[go_p$FDR < alpha]))

# Normalised edata
edata.cpm.df <- edata.cpm %>% 
                add_column(gene_ids = factor(row.names(edata.cpm), 
                                             levels = gene_clusters_lfc), .before = 1) %>%
                gather(key = "sample_name", value = "norm_counts", convert = FALSE, -gene_ids) %>%
                data.frame(., stringsAsFactors = FALSE)

idx <- match(edata.cpm.df$sample_name, design$sample_name)
edata.cpm.df <- cbind(edata.cpm.df, design[idx, !str_detect(colnames(design), "sample_name")])
edata.cpm.df$genotype <- factor(edata.cpm.df$genotype, levels = genotype$name)
#levels(edata.cpm.df$genotype) <- genotype$name
edata.cpm.df$treatment <- factor(edata.cpm.df$treatment[idx], levels = treatment$name)
edata.cpm.df$gene_sig <- as.factor(ifelse(as.character(edata.cpm.df$gene_ids) %in% gene_sig, 1, 0))

logFC_P$genotype <- factor(logFC_P$genotype, levels = genotype$short)
logFC_P$treatment <- factor(logFC_P$treatment, levels = treatment$name)
logFC_P$gene_sig <- as.factor(ifelse(logFC_P$PValue < alpha, 1, 0))

# Plot log2 normalised counts for the GO:TERM --use mclapply for parallel run
mclapply(1:length(go_sig), function(i){
    
    #i = 1
    tmp <- go_sig[i]
    go_gene_sig <- as.character(unique(call_biomart(tmp, 
                                                columns = c("ensembl_gene_id"), 
                                                filter = "go")$ensembl_gene_id))

    nam <- paste(unique(call_biomart(tmp, 
                              columns = c("name_1006", "go_id"), 
                              filter = "go")$name_1006), collapse = ";")

    message(paste("[",i,"/",length(go_sig),"]", sep = ""))
    message(paste("[",tmp,"] Relative normalised counts for differentially expressed genes: ", paste(go_gene_sig, collapse = ";"), " ...", sep = ""))
    
    if(length(go_gene_sig) == 1) {
        message("No. of genes are less")
        flag <- 1
    } else if(length(go_gene_sig) < 1){
        message("NONE GENES FOUND IN ENSEMBL::BIOMART")
        flag <- 2
    } else { 
        flag <- 0
    }
    
    if(flag <= 1){
        
        # CPM Plot
        edata.cpm.df %>% dplyr::filter(gene_ids %in% go_gene_sig) %>%
        ggplot(aes(x = genotype, y = norm_counts, group = treatment)) +
            ggtitle(paste("[",tmp, "] - No. of genes = ", length(go_gene_sig), sep = "")) +
            geom_boxplot(alpha = 0.4, outlier.fill = NA, aes(fill = treatment)) +
            geom_point(position = position_jitter(), 
                       alpha = 0.8,
                       aes(shape = gene_sig, colour = genotype)) +
            facet_grid(~ genotype, scale = "free_x", switch = "x",
                       labeller = labeller(genotype = label.idx)) +
            scale_colour_manual(values = genotype$col.idx, labels = label.idx) +
            scale_fill_manual(values = treatment$col.idx) +
            scale_shape_manual(values = c(`0` = 2, 
                                         `1` = 4), 
                              labels = c(`0` = "Not Significant", 
                                         `1` = "Significant")) +
            theme_AKB +
            theme(axis.text.x = element_blank(),
                  axis.text.y = element_text(size = 16, vjust = 0.5, hjust = 1),
                  axis.title.x = element_text(size = 2, vjust = 1, hjust = 0.5),
                  strip.text.x = element_text(size = 16, 
                                              angle = 60, 
                                              vjust = 0.5, 
                                              hjust = 0.5),
                  legend.text = element_text(size = 12)) +
            labs(x = nam,
                 y = "log2(CPM)",
                 size = "", 
                 fill = "", 
                 colour = "", 
                 shape = "") +
            ggsave(paste(fig.path, "/GSET/GO_CAMERA_genotype/", tmp,"_CPM_sig[",alpha,"]_ngenes[",length(go_gene_sig),"]_genotype.png", sep = ""), 
                   dpi = 300, 
                   device = "png", 
                   height = 16, 
                   width = 9, limitsize = F);

        geno <- which(genotype$short %in% unique(logFC_P$genotype))
                      
        # LFC Plot
        logFC_P %>% dplyr::filter(gene_ids %in% go_gene_sig) %>% 
        arrange(desc(logFC)) %>%
        ggplot(aes(x = treatment, y = logFC)) +
            ggtitle(paste("[",tmp, "] - No. of genes = ", length(go_gene_sig), sep = "")) +
            geom_hline(yintercept = 0, alpha = 1, size = 2) +
            geom_boxplot(alpha = 0.4, outlier.fill = NA, aes(fill = treatment)) +
            geom_point(position = position_jitter(), 
                       alpha = 0.8,
                       aes(shape = gene_sig, colour = genotype)) +
            facet_grid(~ genotype, scale = "free_x", switch = "x",
                       labeller = labeller(genotype = label.idx)) +
            scale_colour_manual(values = genotype$col.idx[geno], labels = label.idx[geno]) +
            scale_fill_manual(values = treatment$col.idx) +
            scale_shape_manual(values = c(`0` = 2, 
                                         `1` = 4), 
                               labels = c(`0` = "Not Significant", 
                                         `1` = "Significant")) +
            theme_AKB +
            theme(axis.text.x = element_text(size = 16, face = "bold", 
                                             angle = 0, vjust = 0.5, hjust = 1),
                  axis.text.y = element_text(size = 16, vjust = 0.5, hjust = 1),
                  strip.text.x = element_text(size = 16, 
                                              angle = 0, 
                                              vjust = 0.5, 
                                              hjust = 0.5, face = "bold"),
                  legend.text = element_text(size = 6), 
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank()) +
            labs(x = "",
                 y = "log2(FC) vs Col-0",
                 size = "", 
                 fill = "", 
                 colour = "", 
                 shape = "") +
            ggsave(paste(fig.path, "/GSET/GO_CAMERA_genotype/", tmp,"_LFC_sig[",alpha,"]_ngenes[",length(go_gene_sig),"]_genotype.png", sep = ""), 
                   dpi = 300, 
                   device = "png", 
                   height = 9, 
                   width = 16, limitsize = F);
        
        # Make data frame for GO separately
        df.tmp <- data.frame(definition = unlist(str_split(nam, ";")))
        df.tmp$go_id <- tmp
        df.tmp$lfc_path <- paste(path, fig.path, "/GSET/GO_CAMERA_genotype/", tmp,"_LFC_sig[",alpha,"]_ngenes[",length(go_gene_sig),"]_genotype.png", sep = "")
        df.tmp$cpm_path <- paste(path, fig.path, "/GSET/GO_CAMERA_genotype/", tmp,"_CPM_sig[",alpha,"]_ngenes[",length(go_gene_sig),"]_genotype.png", sep = "")
        df.tmp$gene_ids <- paste(go_gene_sig, collapse = ";")
        write.table(df.tmp, paste(out.path, "/gset_analysis/GO_CAMERA_genotype/", tmp,"_LFC_sig[",alpha,"]_ngenes[",length(go_gene_sig),"]_genotype.png", sep = ""), sep = "\t")
        
    }else {
        message("Trying next ...")
    }
}, mc.cores = 24)
message("DONE")

# END OF SCRIPT
sessionInfo()
