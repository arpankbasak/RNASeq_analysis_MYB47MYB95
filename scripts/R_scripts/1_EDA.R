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

# Set target for QC
target <- c("AT1G52410", "AT4G39320", 
            "AT5G56870", "AT1G18710", 
            "AT1G74430", "AT5G24780", 
            "AT5G24770", "AT5G24000")

# Raad in data
#edata <- read.table(paste(out.path, "/summarised_exp/edata_rs.txt", sep = ""),
#                       header = TRUE, as.is = TRUE)
load(paste(out.path, "/summarised_exp/salmon_summarized_experiments.Rdata", sep = ""))
load(paste(out.path, "/summarised_exp/salmon_summarized_experiments_genelevel.Rdata", sep = ""))
edata <- as.data.frame(tximport_sum$counts)
edata.length <- as.data.frame(tximport_sum$length)
# edata <- as.data.frame(tximport_obj$counts)
# edata.length <- as.data.frame(tximport_obj$length)


metadata <- read.table(paste(data.files, "metadata.txt", sep = ""),
                       header = TRUE,   as.is = TRUE)

metadata$genotype <- str_replace(metadata$genotype, "-tKO", "2/3/4")
metadata$genotype_short <- genotype$short[match(metadata$genotype, genotype$name)]
#metadata$sample_name <- str_replace(metadata$sample_name, "_WT", "WT")
#colnames(edata) <- str_replace(colnames(edata), pattern = ".genome.bam", replacement = "")

levels_group <- sapply(genotype[, "short"], function(x) paste(x, treatment[, "name"], sep = "_"))

design <- metadata %>% 
    group_by(sample_name, treatment, genotype, genotype_short, random, replicate, seq_run) %>% 
    summarise(temp = n()) %>%
    dplyr::select(-temp) %>%
    mutate(group = paste(genotype_short, treatment, sep = "_")) %>%
    data.frame(.)

# design <- design %>% filter(!sample_name %in% c("myb47_95_DW5_27")) %>% data.frame(.)

# design$treatement <- factor(design$treatment, levels = treatment[, "name"])
# design$genotype <- factor(design$genotype, levels = genotype[, "name"], labels = genotype[, "short"])
design$group <- factor(design$group, levels = levels_group[levels_group %in% design$group])
# design$random <- as.factor(paste("B", design$random, design$random, sep = ""))
# design$replicate <- as.factor(design$replicate)
design$seq_run <- as.factor(design$seq_run)

write.table(design,
            file = paste(out.path, "/design.txt", sep = ""),
            sep = "\t")

# Clean edata columns
# idx <- which(colnames(edata) %in% design$sample_name)
# edata <- edata[,idx]

# idx <- which(colnames(edata.length) %in% design$sample_name)
# edata.length<- edata.length[,idx]
# # If sample to be replaced ----------------------------------------------------
# ix <- which(design$sample_name == "myb47_95_DW5_27")
# iy <- which(design$sample_name == "myb47_95JA6_31")
# design$sample_name[c(ix, iy)] <- c("myb47_95JA6_31", "myb47_95_DW5_27")
# # ------------------------------------------------------------------------------

# Align NP matrix and sample groups
row.names(design) <- design$sample_name
design <- design[colnames(tximport_sum$counts), c("group", "seq_run")] 

# Remove the effect of sequencing run
design$seq_run <- droplevels(design$seq_run)
design$random <- droplevels(design$random)
design.matrix <- model.matrix(~ 1 - group + seq_run, data = design)
null.matrix <- model.matrix(~ 1 - seq_run, data = design)
# design.matrix <- model.matrix(~ 0 + seq_run, data = design)
colnames(design.matrix) <- str_replace_all(colnames(design.matrix), "group|random", "")
design.matrix <- design.matrix[, colSums(design.matrix) != 0]
design.matrix <- design.matrix[rowSums(design.matrix) != 0,]
                       
# Normalisation by RLOG method
# dds <- DESeqDataSetFromTximport(tximport_sum, design, ~seq_run + group)
dds <- DESeqDataSetFromTximport(tximport_sum, design, design.matrix)
# If sample to be removed ----------------------------------------------------
# ex <- c("myb47_95_DW5_27", "myb47_95JA6_31", 
#         "WT_JA_4_28", "WT_DW_3",
#         "mycTKO_JA_13", "mycTKO_DW_6_20")

# ix <- which(!colData(dds)$sample_name %in% ex)
# design <- as.data.frame(colData(dds)[ix,])
# ix <- which(!colnames(assay(dds)) %in% ex)
# dds <- dds[,ix]
# # ------------------------------------------------------------------------------

# Remove the effect of sequencing run
# design$seq_run <- droplevels(design$seq_run)
# design$random <- droplevels(design$random)
# design.matrix <- model.matrix(~ 0 + group + seq_run, data = design)
# null.matrix <- model.matrix(~ 0 + seq_run, data = design)
# # design.matrix <- model.matrix(~ 0 + seq_run, data = design)
# colnames(design.matrix) <- str_replace_all(colnames(design.matrix), "group|random", "")
# design.matrix <- design.matrix[, colSums(design.matrix) != 0]
# design.matrix <- design.matrix[rowSums(design.matrix) != 0,]

# design(dds) <- design.matrix
contrast.mat <- limma::makeContrasts(
     
    # Effect of treatment
    treatment_WT_JA = (WT_JA - WT_DW),
    treatment_myb_JA = (myb_JA - myb_DW),
    treatment_myc.tKO_JA = (myc.tKO_JA - myc.tKO_DW),
    
    # Effect of genotype
    genotype_WT_myb_JA = (myb_JA - WT_JA),
    genotype_WT_myb_DW = (myb_DW - WT_DW),
    genotype_WT_myc.tKO_JA = (myc.tKO_JA - WT_JA),
    genotype_WT_myc.tKO_DW = (myc.tKO_DW - WT_DW),
    
    levels = design.matrix)

contrast.names <- attr(contrast.mat, "dimnames")$Contrasts
n <- length(contrast.names)

# Filter lowcounts
keep <- rowSums(counts(dds) >= 10) >= 1
dds <- dds[keep,]
                       
#dds <- DESeqDataSetFromMatrix(as.matrix(edata), 
#                             colData = design,
#                             design = ~ 0 + genotype + treatment)
# Run DESEQ pipeline
dds <- DESeq2::estimateSizeFactors(dds)
dds <- DESeq(dds, test = "LRT", reduced=null.matrix, betaPrior=FALSE)
dds_all <- dds
# dds <- DESeq(dds)

# combi_group <- colnames(rowData(dds))[str_detect(colnames(rowData(dds)), "^group_")]

dds.gene.list <- lapply(contrast.names, function(x) {
    
    dds1 <- results(dds_all, 
            contrast= list("group", 
              colnames(contrast.mat)[contrast.mat[, which(contrast.names == i)] == 1], 
              colnames(contrast.mat)[contrast.mat[, which(contrast.names == x)] == -1]), 
            pAdjustMethod="BH", 
#             lfcThreshold=FC_threshold, 
            test="LRT", 
            cooksCutoff = TRUE,
            independentFiltering= FALSE,
#             addMLE = TRUE,
            alpha = alpha)
    
    results(dds_all, 
            contrast= list("WT_JA", "myb_DW"), 
            pAdjustMethod="BH", 
#             lfcThreshold=FC_threshold, 
            test="LRT", 
            cooksCutoff = TRUE,
            independentFiltering= FALSE,
#             addMLE = TRUE,
            alpha = alpha)

#   dds1 <- results(dds_all, 
#             contrast= contrast.mat[, which(contrast.names == x)], 
#             pAdjustMethod="BH", 
# #             lfcThreshold=FC_threshold, 
#             test="LRT", 
#             cooksCutoff = TRUE,
#             independentFiltering= FALSE,
# #             addMLE = TRUE,
#             alpha = alpha)

    return(dds1)
}) # lfcThreshold="if required", filter="lowcount", test="LRT"

temp <- as.data.frame(rowData(dds))
dds.gene.list <- mclapply(row.names(temp), function(x){

  temp.gene <- temp[x,]
  treatment_WT_JA <- binom.test(temp.gene$WT_JA, temp.gene$WT_DW)

})
names(dds.gene.list) <- contrast.names

# LogFC and PValue tables
logFC_P.list <- lapply(1:n, function(x) {
  table <- data.frame(logFC = dds.gene.list[[x]][,"log2FoldChange"],
                     PValue =  p.adjust(dds.gene.list[[x]][,"padj"]))
  #table <- table %>% dplyr::mutate(PValue = p.adjust(table[,2], method=p.adj.method))
  colnames(table) <- paste(contrast.names[x], colnames(table), sep="_")
    row.names(table) <- row.names(dds.gene.list[[x]])
  return(table)
})

# filtered <- lapply(1:n, function(x) { 
#     idx <- which(is.na(dds.gene.list[[x]][,"padj"]) & 
#                       dds.gene.list[[x]][,"baseMean"] > 0)
#     idx <- row.names(dds.gene.list[[x]])[idx]
# return(idx)
# })
# names(filtered) <- contrast.names
# filtered <- do.call(data.frame, filtered)
                       
# Make LFC Table
logFC_P <- do.call(data.frame, logFC_P.list)
logFC_P.cook <- cbind.data.frame(logFC_P, rowData(dds)[row.names(logFC_P),c("dispOutlier", "maxCooks")])
logFC_P.long <- logFC_P %>%
                       add_column(gene_id = row.names(logFC_P), .before = 1) %>%
                       gather(key = "key", value = "val", convert = FALSE, -gene_id) %>%
                       data.frame(.)
                       
# row.names(logFC_P) <- gene.ids
DE.list <- lapply(names(dds.gene.list), function(x){
    table <- data.frame(logFC = dds.gene.list[[x]][,"log2FoldChange"],
                        PValue = dds.gene.list[[x]][,"padj"]) %>% 
    mutate(decided = ifelse(PValue < alpha & abs(logFC) > FC_threshold, 1 * sign(logFC), 0)) %>% 
    dplyr::select(decided)
    colnames(table) <- paste(x)
    row.names(table) <- row.names(dds.gene.list[[x]])
  return(table)
})

names(DE.list) <- names(dds.gene.list)           
DEGs <- do.call(data.frame, DE.list)
DEGs <- as.data.frame(apply(DEGs, 2, function(x) replace_na(x, 0)))
count.mat <- data.frame(total = sapply(DEGs, function(x) sum(abs(x))),
                        upregulated = sapply(DEGs, function(x) sum(x ==  1)),
                        downregulated = sapply(DEGs, function(x) sum(x == -1)))
                                               
# rld <- rlog(dds, blind = FALSE)
rld <- rlog(dds)
log.norm <- as.data.frame(assay(rld))
z.norm <- as.data.frame(apply(assay(dds), 2, function(x) (x - mean(x))/sd(x)))
row.names(z.norm) <- row.names(as.data.frame(assay(dds)))
mean.norm <- as.data.frame(rowData(dds)[, as.character(design$group)])

                       
# QC fro library size
# df <- data.frame(x = sizeFactors(dds), 
#                  y = colSums(counts(dds)),
#            genotype = design$genotype[match(names(sizeFactors(dds)), design$sample_name)],
#            treatment = design$treatment[match(names(sizeFactors(dds)), design$sample_name)])

# df$genotype <- factor(df$genotype, levels = genotype$name)
# df$treatment <- factor(df$treatment, levels = treatment$name)

# message("Spilling size factor estimation plot by rLog method ...")
# df %>% ggplot(aes(x, y)) +
#     geom_point(alpha = 0.6, aes(colour = treatment)) +
#     geom_smooth(method = "lm", size = 0.5, colour = "red") +
#     scale_colour_manual(values = treatment$col.idx) +
#     facet_grid(.~ genotype, scales = "free", labeller = labeller(genotype = label.idx)) +
#     theme_AKB +
#     labs(x = "", y = "", colour = "") + 
#     ggsave(paste(fig.path, "/EDA_sizefactor_estimation.png", sep = ""), 
#            dpi = 300, device = "png", 
#            height = 10, 
#            width = 10, limitsize = F)

message("DONE")

write.table(log.norm,
            file = paste(out.path, "/edata.log.normalised_DESeq2.txt", sep = ""),
            sep = "\t")

save(list = c("design", "dds", "rld", "logFC_P", "logFC_P.cook", "logFC_P.long", "DEGs", "count.mat"), 
  file = paste(out.path, "/DESeq2_objects.Rdata", sep = ""))
write.table(logFC_P, paste(stat.path, "/logFC_P.deseq.txt", sep = ""), sep = "\t")
write.table(logFC_P.long, paste(stat.path, "/logFC_P.long.deseq.txt", sep = ""), sep = "\t")
write.table(logFC_P.cook, paste(stat.path, "/logFC_P_with_oitlier_genes.deseq.txt", sep = ""), sep = "\t")
write.table(DEGs,
            paste(stat.path, "/DEGs.deseq.txt", sep = ""),
            sep = "\t")
write.table(count.mat,
            paste(stat.path, "/up_down_reg.deseq.txt", sep = ""),
            sep = "\t")

# Make clusters of samples and genes -- this will be used as factor levels
d <- 1-(cor(log.norm))
clust_sample <- hclust(as.dist(d), method = "ward.D")
sample_clusters <- colnames(log.norm)[clust_sample$order] # Cluster grouping is required by concensus clustering method
mds_obj <- cmdscale(d, k = 3, eig = TRUE)

d <- 1 - cor(t(log.norm))
clust_gene <- hclust(as.dist(d), method = "ward.D")
gene_clusters <- row.names(log.norm)[clust_gene$order]
rm(list = c("d")) # save RAM

# MDS Plot considerations
eig <- mds_obj$eig/sum(mds_obj$eig)
idx <- match(design$sample_name, row.names(mds_obj$points))
df <- cbind(design, 
            mds1 = mds_obj$points[idx,1], 
            mds2 = mds_obj$points[idx,2], 
            mds3 = mds_obj$points[idx,3])

message("Spilling MDS Plot ...")
df %>% mutate(sample_name = factor(sample_name, levels = sample_clusters)) %>%
        ggplot(aes(x = saturate(mds1), y = saturate(mds2))) +
        geom_hline(yintercept = 0.0, linetype = "solid", colour = "darkgrey", alpha = 0.6, size = 0.8) +
        geom_vline(xintercept = 0.0, linetype = "solid", colour = "darkgrey", alpha = 0.6, size = 0.8) +
        geom_point(alpha = 0.6, aes(shape = treatment, colour = genotype), size = 5) +
#         geom_text(aes(colour = genotype, label = sample_name), 
#                   nudge_x = 0.5, nudge_y = 0.5,
#                   check_overlap = TRUE, 
#                   show.legend = FALSE, 
#                   size = 1) +
        scale_colour_manual(values = genotype$col.idx, labels = label.idx) +
        scale_shape_manual(values = treatment$pch.idx) +
        stat_ellipse(level = 0.95, lwd = 1, type = "t", aes(group = treatment, linetype = treatment)) +
        scale_linetype_manual(values = treatment$lty.idx) +
        # xlim(c(-1,1)) +
        # ylim(c(-1,1)) +
        theme_AKB +
        labs(x = paste("MDS 1: ", round(eig[1], 4) * 100, " %", sep = ""), 
             y = paste("MDS 2: ", round(eig[2], 4) * 100, " %", sep = ""), 
             colour = "",
             shape = "") +
        theme(legend.title = element_blank(), 
              legend.text = element_text(size = 6),
              axis.title = element_text(hjust = 0.5, vjust = 0.5),
              panel.border = element_rect(size = 1, 
                fill = "transparent"), 
              axis.line = element_blank()) +
        ggsave(paste(fig.path, "/EDA_MDS_Plot.png", sep = ""), 
               dpi = 600, 
               device = "png", 
               units = "in",
               height = 5, 
               width = 5, limitsize = F, bg = "transparent")
message("DONE")

# Plot Canvas

message("Spilling Relative expression heatmap ...")
# Plot canvas for the relative expression profile for all
df <- as.data.frame(z.norm) %>% 
    add_column(gene_ids = row.names(z.norm), .before = 1) %>%
    gather(key = "sample_name", value = "counts", convert = FALSE, -gene_ids) %>%
    mutate(sample_name = factor(sample_name, levels = sample_clusters),
          gene_ids = factor(gene_ids, levels = gene_clusters)) # Add a pseudocounts

df$genotype <- design$genotype[match(df$sample_name, design$sample_name)]
df$treatment <- design$treatment[match(df$sample_name, design$sample_name)]

df$genotype <- factor(df$genotype, levels = genotype$name)
df$treatment <- factor(df$treatment, levels = treatment$name)

levels(df$genotype) <- genotype$short

message("Conducting Cannonical Correlation analysis supported by permutation test ...")
levels <- unlist(lapply(levels(df$treatment), function(x) paste(levels(df$genotype), x, sep = "_")), use.names = FALSE)

df <- df %>% mutate(group = factor(paste(genotype, treatment, sep = "_"), 
                                       levels = levels), 
                    random = design$seq_run[match(as.character(sample_name), row.names(design))])


# Fetch the factor of interest from the data
f <- formula(t(log.norm) ~ group + Condition(seq_run))
cca_obj <- vegan::cca(formula = f, 
                    design, distance = "pearson", sqrt.dist = T, add = F)

# Conduct PERMANOVA
set.seed(1)
p.nova <- vegan::anova.cca(cca_obj)

# Extract variables from stat objects
eig <- format(100 * (cca_obj$CCA$eig/sum(cca_obj$CCA$eig)), digits = 4)
pval <- p.nova$`Pr(>F)`[1]
chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
                       row.names = c("total", "constrianed", "unconstrained"))

# Format P-Value for representation
variance <- format(100 * variable$proportion[2], digits = 4)
ti <- paste("Genotype + Treatment = ", variance, "% p = ", pval, 
             ifelse(pval < 0.05, "*", " (NS)"), sep = "")

message("Spilling Canvas for CCA analysis ...")
# Plot canvas

design %>% 
    dplyr::mutate(replicate = as.factor(replicate), 
                  CCA1 = as.numeric(cca_obj$CCA$wa[row.names(design),1]),
                  CCA2 = as.numeric(cca_obj$CCA$wa[row.names(design),2])) %>% 
    ggplot2::ggplot(aes(x= saturate(CCA1), y = saturate(CCA2))) +
    ggtitle(ti) + 
    geom_point(alpha = 0.6, 
               size = 5, 
               aes(colour = genotype, shape = treatment)) +
    geom_hline(yintercept = 0.0, linetype = "dashed", colour = "darkgrey", alpha = 0.6, size = 0.8) +
    geom_vline(xintercept = 0.0, linetype = "dashed", colour = "darkgrey", alpha = 0.6, size = 0.8) +
    stat_ellipse(type = "t", alpha = 0.6, linetype = "dashed",
                 level = 0.95, aes(group = genotype, 
                                   colour = genotype)) +
    scale_color_manual(values = genotype$col.idx, labels = label.idx) + 
    scale_shape_manual(values = treatment$pch.idx) +
    # scale_linetype_manual(values = genotype$lty.idx, guide = FALSE) +
    theme_AKB +
    theme(legend.title = element_blank(), 
              legend.text = element_text(size = 6),
              axis.title = element_text(hjust = 0.5, vjust = 0.5),
              panel.border = element_rect(size = 1, 
                fill = "transparent"), 
              axis.line = element_blank()) +
    labs(x = paste("CCA 1: ",eig[1],"%", sep = ""), 
         y = paste("CCA 2: ",eig[2],"%", sep = ""), 
         colour = "", shape = "", lty = "") +
    ggsave(filename = paste(fig.path,"/cca_analysis.png", sep = ""), 
           dpi = 600, 
               device = "png", 
               units = "in",
               height = 5, 
               width = 5, limitsize = F, bg = "transparent")


# # Statistics for genewise analysis on the target
# res <- mclapply(1:length(target), function(x){

#     # x <- 1
#     # Take out gene
#     i <- target[x]

#     # Slice out from the matrix
#     dat <- log.norm[i, ] %>% gather(key = "sample_id", value = "log.norm", convert = FALSE) %>%
#     mutate(genotype = design$genotype[match(sample_id, design$sample_name)],
#            treatment = design$treatment[match(sample_id, design$sample_name)],
#            seq_run = as.factor(paste0("B", design$seq_run[match(sample_id, design$sample_name)]))) %>%
#     mutate(group  = as.factor(paste(genotype, treatment, sep = "_"))) %>%
#     data.frame(., stringsAsFactors = FALSE)


#     # Do statistical testing using generalised linear model
#     mod <- aov(glm(log.norm ~ 0 + group + seq_run, data = dat))
#     mod <- TukeyHSD(mod)
#     df <- as.data.frame(mod$group) %>% 
#     add_column(geneid = i, x = row.names(.), .before = 1) %>%
#     separate(x, into = c("group1", "group2"), 
#         sep = "_", convert = F, remove = F) %>%
#     mutate(group1 = as.factor(group1), 
#         group2 = as.factor(group2), 
#         significance = ifelse(`p adj` < 0.05, "*", "")
#         )

#     return(df)

# }, mc.cores = 6)

# res.df <- do.call(rbind.data.frame, res)

# res.df %>% write.table(., paste(stat.path, "/tsd_statistics_target.txt", sep = ""), sep = "\t")

# mclapply(1:length(target), function(x){

#     # x <- 8
#     # Take out gene
#     i <- target[x]

#     # Slice out from the matrix
#     dat <- log.norm[i, ] %>% gather(key = "sample_id", value = "log.norm", convert = FALSE) %>%
#     mutate(genotype = design$genotype[match(sample_id, design$sample_name)],
#            treatment = design$treatment[match(sample_id, design$sample_name)],
#            seq_run = as.factor(paste0("B", design$seq_run[match(sample_id, design$sample_name)]))) %>%
#     mutate(group  = as.factor(paste(genotype, treatment, sep = "_"))) %>%
#     data.frame(., stringsAsFactors = FALSE)

#     # Do statistical testing using generalised linear model
#     dat %>% ggplot(aes(y = log.norm, x = genotype)) +
#     ggtitle(i) +
#     geom_point(aes(colour = treatment, shape = genotype), position = "jitter", alpha = 0.5) +
#     geom_boxplot(aes(colour = treatment), fill = "transparent", alpha = 0.8) +
#     scale_colour_manual(values = treatment$col.idx) +
#     scale_shape_manual(values = genotype$pch.idx, labels = label.idx) +
#     theme_AKB +
#     theme(axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5, size = 10),
#           axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
#           strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
#           legend.text = element_text(size = 6)) +
#     labs(x = "", colour = "", shape = "", y = "log.normalised[counts]") +
#     ggsave(paste(fig.path, "/hallmark_", i,".png", sep = ""),
#            dpi = 300,
#            device = "png",
#            height = 9,
#            width = 16, limitsize = F)
    
#     dat %>% ggplot(aes(y = log.norm, x = sample_id)) +
#     ggtitle(i) +
#     geom_point(aes(colour = treatment, shape = genotype)) +
# #     geom_boxplot(aes(colour = treatment), fill = "transparent", alpha = 0.8) +
#     scale_colour_manual(values = treatment$col.idx) +
#     scale_shape_manual(values = genotype$pch.idx, labels = label.idx) +
#     theme_AKB +
#     theme(axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5, size = 10),
#           axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
#           strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
#           legend.text = element_text(size = 6)) +
#     labs(x = "", colour = "", shape = "", y = "log.normalised[counts]") +
#     ggsave(paste(fig.path, "/hallmark_test", i,".png", sep = ""),
#            dpi = 300,
#            device = "png",
#            height = 9,
#            width = 16, limitsize = F)

# })


message("DONE")


# Save Objects
save(list = c("clust_sample", "clust_gene"), file = paste(out.path, "/clusters.Rdata", sep = ""))
save(list = "df", file = paste(out.path, "/log.normalised.matrix.Rdata", sep = ""))

# END OF SCRIPT
sessionInfo()
