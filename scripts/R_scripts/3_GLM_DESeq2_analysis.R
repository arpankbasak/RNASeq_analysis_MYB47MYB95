rm(list = ls())
# NOTE::
# USE DESEQ2 module
options(warn = 1,
       mc.cores = 8)

pkgs <- c("tidyverse", "biocparallel", "DESeq2", "limma", "edgeR", "csaw")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)
# Raad in data
load(paste(out.path, "/summarised_exp/salmon_summarized_experiments_genelevel.Rdata", sep = ""))
edata <- as.data.frame(tximport_sum$counts)
edata.length <- as.data.frame(tximport_sum$length)
#colnames(edata) <- str_replace(colnames(edata), pattern = ".genome.bam", replacement = "")

# TRy to remove low conts
idx <- cpm(edata)
idx <- rowSums(idx) >= 10
edata.filtered <- edata[idx,]
edata.length.filtered <- edata[row.names(edata.filtered),]

write.table(edata.filtered, 
            paste(out.path, "/summarised_exp/edata.filtered_rs.txt", sep = ""), 
            sep = "\t")

design <- read.table(paste(out.path, "/design.txt", sep = ""),
                       header = TRUE, as.is = TRUE)
# row.names(design) <- design$sample_name

load(paste(out.path, "/DESeq2_objects.Rdata", sep = ""))
# design <- design %>% mutate(random = paste("B", random, sep = ""), 
                            # seq_run = paste("B", seq_run, sep = ""))
# design$genotype <- genotype$short[match(design$genotype, genotype$name)]
# design <- design %>% mutate(group = paste(genotype, treatment, sep = "_"))

# levels <- c(sapply(genotype$short, paste, treatment$name, sep = "_"))
groups <- factor(design$group, levels = levels(design$group))
random <- as.factor(design$seq_run)

model <- model.matrix(~ 0 + groups + random) # Drop intercept
colnames(model) <- str_replace(colnames(model), "groups", "")
model <- model[,which(colSums(model) != 0)]

edata.filtered <- edata.filtered[,design$sample_name]

# Build contrast
contrast.mat <- makeContrasts(
     
    # Effect of treatment
    treatment_WT_JA = (WT_JA - WT_DW),
    treatment_myb_JA = (myb_JA - myb_DW),
    treatment_myc.tKO_JA = (myc.tKO_JA - myc.tKO_DW),
    
    # Effect of genotype
    genotype_WT_myb_JA = (myb_JA - WT_JA),
    genotype_WT_myb_DW = (myb_DW - WT_DW),
    genotype_WT_myc.tKO_JA = (myc.tKO_JA - WT_JA),
    genotype_WT_myc.tKO_DW = (myc.tKO_DW - WT_DW),
    
    # Effect of genotype normalised to treatment
    genotypetreatment_WT_myb_JA = (myb_JA - myb_DW) - (WT_JA - WT_DW),
    genotypetreatment_WT_myc.tKO_JA = (myc.tKO_JA - myc.tKO_DW) - (WT_JA - WT_DW),

    
    levels = model)

# Store contrast names for later use
contrast.names <- attr(contrast.mat, "dimnames")$Contrasts
n <- length(contrast.names)

# Prior normalisation per length
edata.norm <- edata.filtered/exp(rowMeans(log(edata.filtered)))
counts.norm <- edata.filtered/edata.length.filtered

# Estimate the effective library size
eff.lib <- calcNormFactors(counts.norm) * colSums(edata.norm)
edata.norm <- sweep(edata.norm, 2, eff.lib, "*")
edata.norm <- log(edata.norm)

delist <- DGEList(counts = edata.filtered, group = groups)
delist <- scaleOffset(delist, as.matrix(edata.norm))
idx <- filterByExpr(delist)
delist <- delist[idx,]

# delist <- DGEList(counts = edata.filtered, group = groups)
# delist <- calcNormFactors(delist)

# Estimate common and tag-wise disperse
de <- estimateGLMCommonDisp(delist, model)
de <- estimateGLMTagwiseDisp(de, model)
gene.ids <- row.names(de$counts)

save(list = c("gene.ids", "de", "delist", 
              "contrast.mat", "contrast.names", 
              "model", "n", "groups"), 
     file = paste(out.path, "edger_objects.GLM.Rdata", sep = ""))

fit <- glmFit(de, model)
# LRT.df <- glmLRT(fit)
# LRT for each contrasts
LRT.list <- lapply(contrast.names, function(x) glmLRT(fit, contrast=contrast.mat[, which(contrast.names == x)]))
names(LRT.list) <- contrast.names

# LogFC and PValue tables
logFC_P.list <- lapply(1:n, function(x) {
  table <- LRT.list[[x]]$table[,c(1,4)]
  table <- table %>% dplyr::mutate(PValue = p.adjust(table[,2], method=p.adj.method))
  colnames(table) <- paste(contrast.names[x], colnames(table), sep="_")
  return(table)
})

logFC_P <- do.call(data.frame, logFC_P.list)
row.names(logFC_P) <- gene.ids
                   
write.table(logFC_P, paste(stat.path, "/logFC_P.edger.txt", sep = ""), sep = "\t")

logFC_P.list <- lapply(1:n, function(x) {
    table <- LRT.list[[x]]$table[,c(1,4)]
    table$PValue <- p.adjust(table[,2], method=p.adj.method)
    table$gene_ids <- row.names(table)
    table$contrast <- contrast.names[x]
  #colnames(table) <- paste(contrast.names[x], colnames(table), sep="_")
  return(table)
})
logFC_P.long <- do.call(rbind, logFC_P.list)
write.table(logFC_P.long, paste(stat.path, "/logFC_P.long.edger.txt", sep = ""), sep = "\t")

# Significance picking for each tested model
contrasts_idx <- intersect(names(LRT.list), contrast.names)
DE.list <- lapply(contrasts_idx, function(x) decideTestsDGE(LRT.list[[x]], 
                                                              adjust.method=p.adj.method, 
                                                              p.value=alpha))
names(DE.list) <- contrasts_idx           
DEGs <- as.data.frame(DE.list)
colnames(DEGs) <- contrasts_idx
row.names(DEGs) <- gene.ids
write.table(DEGs,
            paste(stat.path, "/DEGs.edger.txt", sep = ""),
            sep = "\t")

# Number of significant differentially expressed genes at 5% FDR
count.mat <- data.frame(total = sapply(DEGs, function(x) sum(abs(x))),
                        upregulated = sapply(DEGs, function(x) sum(x ==  1)),
                        downregulated = sapply(DEGs, function(x) sum(x == -1)))

write.table(count.mat,
            paste(stat.path, "/up_down_reg.edger.txt", sep = ""),
            sep = "\t")
                                                                           
log2cpm <- cpm(de, 
               prior.count = 2, 
               log = T) # Counts per million in log2 scale

write.table(log2cpm, paste(stat.path, "/log2cpm.edger.txt", sep = ""), sep = "\t")

                                        
                   
save(list = c("logFC_P", "logFC_P.long", "DEGs", "log2cpm", "design"), 
     file = paste(out.path, "edger_objects.plotting.Rdata", sep = ""))

# END OF SCRIPT
sessionInfo()