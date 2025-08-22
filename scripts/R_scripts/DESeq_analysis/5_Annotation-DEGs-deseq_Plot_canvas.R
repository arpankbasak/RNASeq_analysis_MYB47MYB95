rm(list = ls())

options(warn = 1,
       mc.cores = 6)

pkgs <- c("tidyverse", "GO.db", "parallel", "org.At.tair.db", "AnnotationDbi", "grid", "gridExtra")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
load(paste(out.path, "/gset_analysis/objects/deseq-lfc_annotated.Rdata", sep = ""))
load(paste(out.path, "/LFC-cluster_deseq-treatment.Rdata", sep = ""))

message("Aracyc and KEGG pathway mapping ...")
idx <- str_detect(colnames(logFC_P.annotated), "^treatment")
df <- logFC_P.annotated[,idx] %>%
   add_column(gene_ids = row.names(.), .before = 1) %>%
   gather(key  = "key", value = "vals", convert = FALSE, -gene_ids) %>%
   separate(key, into = c("x1", "Genotype", "x2", "idd"), sep = "_") %>%
   spread(key = idd, value = vals, convert = FALSE, fill = NA) %>%
   mutate(Genotype = factor(Genotype, levels = genotype$short),
      aracyc = logFC_P.annotated$aracyc[match(.$gene_ids, row.names(logFC_P.annotated))],
      kegg_path = logFC_P.annotated$kegg_path[match(.$gene_ids, row.names(logFC_P.annotated))],
      gene_name = logFC_P.annotated$gene_name[match(.$gene_ids, row.names(logFC_P.annotated))],
      sig = ifelse(PValue < alpha & abs(logFC) > FC_threshold, 1, 0)
      )

# Map ARACYC
# df$term_aracyc <- ""
# df$term_aracyc[which(!str_detect(df$aracyc, "NA"))] <- "Annotated"

# # Kegg Pathway annotated
# df$term_kegg <- ""
# df$term_kegg[which(!str_detect(df$kegg_path, "NA"))] <- "Annotated"

# Consider ARACYC terms
# df[,c("gene_ids", "term_aracyc")] %>%
#          mutate(gene_ids = factor(as.character(gene_ids), levels = gene_clusters_lfc),
#             tag = as.factor(term_aracyc), 
#             temp = as.factor("x")
#             ) %>%
#           ggplot(aes(x = temp, 
#                         y = gene_ids, 
#                         fill = tag), na.rm = FALSE) +
#           geom_tile(alpha = 1) +
#           scale_fill_manual(values = c(`Annotated` = "red", 
#             `''` = "transparent"), guide = FALSE) +
#           theme_AKB +
#           theme(panel.spacing = unit(1, "lines"),
#              panel.border = element_rect(size = 0.3, fill = "transparent"),
#                axis.text.x = element_blank(),
#                 axis.text.y = element_blank(),
#                 axis.line.x = element_blank(),
#                 axis.line.y = element_blank(),
#                 axis.ticks= element_blank(),
#                 legend.text = element_text(size = 6)) +
#           labs(x = "", y = "",
#               fill = "") +
#           ggsave(filename = paste(fig.path, "/DEGs/maps_aracyc.png", sep = ""), 
#            dpi = 600, 
#            device = "png", 
#            height = 21, 
#            width = 3, limitsize = F, bg = "transparent")


# Consider KEGG terms
# df[,c("gene_ids", "term_kegg")] %>%
#          mutate(gene_ids = factor(as.character(gene_ids), levels = gene_clusters_lfc),
#             tag = as.factor(term_kegg), 
#             temp = as.factor("x")
#             ) %>%
#           ggplot(aes(x = temp, 
#                         y = gene_ids, 
#                         fill = tag), na.rm = FALSE) +
#           geom_tile(alpha = 1) +
#           scale_fill_manual(values = c(`Annotated` = "red", 
#             `''` = "transparent"), guide = FALSE) +
#           theme_AKB +
#           theme(panel.spacing = unit(1, "lines"),
#              panel.border = element_rect(size = 0.3, fill = "transparent"),
#                axis.text.x = element_blank(),
#                 axis.text.y = element_blank(),
#                 axis.line.x = element_blank(),
#                 axis.line.y = element_blank(),
#                 axis.ticks= element_blank(),
#                 legend.text = element_text(size = 6)) +
#           labs(x = "", y = "",
#               fill = "") +
#   ggsave(filename = paste(fig.path, "/DEGs/maps_kegg.png", sep = ""), 
#            dpi = 600, 
#            device = "png", 
#            height = 21, 
#            width = 3, limitsize = F, bg = "transparent")

# Mark the JA responsive genes
# df$term_response <- ""
# df$term_response[str_detect(df$gene_name, "insect|herbivore")] <- "Herbivory"
# df$term_response[str_detect(df$gene_name, "jasmonate|jasmonic acid")] <- "JA response"
# df$term_response[str_detect(df$gene_name, "glucosinolate")] <- "Glucosinolate"
# df$term_response[str_detect(df$gene_name, "ER body")] <- "ER Body"
# df$term_response[str_detect(df$gene_name, "insect|herbivore") & str_detect(df$gene_name, "jasmonate")] <- "Herbivory & JA response"
# df$term_response[str_detect(df$gene_name, "glucosinolate") & str_detect(df$gene_name, "insect|herbivore")] <- "Herbivory & Glucosinolate"
# df$term_response[str_detect(df$gene_name, "glucosinolate") & str_detect(df$gene_name, "insect|herbivore")] <- "JA response & Glucosinolate"

# Annotate the JA responsive genes
# df[,c("gene_ids", "term_response")] %>%
#          mutate(gene_ids = factor(as.character(gene_ids), levels = gene_clusters_lfc),
#             tag = as.factor(term_response), 
#             temp = as.factor("x")
#             ) %>%
#           ggplot(aes(x = temp, 
#                         y = gene_ids, 
#                         fill = tag), na.rm = FALSE) +
#           geom_tile(alpha = 1) +
#           scale_fill_manual(values = c(
#             `Herbivory` = "darkolivegreen", 
#             `JA response` = "darksalmon", 
#             `Herbivory & JA response` = "darkmagenta", 
#             `ER body` = "darkturquoise", 
#             `Glucosinolate` = "darkorange", 
#             `Herbivory & Glucosinolate` = "darkviolet", 
#             `JA response & Glucosinolate` = "darkorchid", 
#             `''` = "transparent")
#           ) +
#           theme_AKB +
#           theme(panel.spacing = unit(1, "lines"),
#             panel.border = element_rect(size = 0.3, fill = "transparent"),
#                axis.text.x = element_blank(),
#                 axis.text.y = element_blank(),
#                 axis.line.x = element_blank(),
#                 axis.line.y = element_blank(),
#                 axis.ticks= element_blank(),
#                 legend.text = element_text(size = 6)) +
#           labs(x = "", y = "",
#               fill = "") +
#   ggsave(filename = paste(fig.path, "/DEGs/maps_genename.png", sep = ""), 
#            dpi = 600, 
#            device = "png", 
#            height = 21, 
#            width = 3, limitsize = F, bg = "transparent")

# Aracyc cumsum PValue Heatmap
# df_aracyc <- df %>% group_by(aracyc, Genotype) %>%
# summarise(pvals = mean(PValue), 
#    sig = sum(PValue < alpha),
#    ngenes = n()) %>%
# filter(ngenes > 2, !is.na(pvals), sig != 0) %>% 
# na.omit(.) %>%
# mutate(sig = sig/ngenes) %>%
# mutate(pvals = -log10(pvals*sig)) %>%
# arrange(pvals, desc(sig)) %>%
# data.frame(.)

# df_aracyc %>%
#          mutate(aracyc = as.factor(aracyc), 
#             ) %>%
#           ggplot(aes(x = pvals, 
#                         y = aracyc), na.rm = FALSE) +
#           geom_point(alpha = 1,
#             aes(colour = pvals, size = ngenes)) +
#           facet_grid(.~ Genotype, switch = "x" ,scale = "free", space = "free") +
#           scale_colour_gradient(low = "blue", high= "darkred") +
#           theme_AKB +
#           theme(panel.spacing = unit(1, "lines"),
#                axis.text.x = element_text(size = 6),
#                 axis.text.y = element_text(size = 4, hjust = 0, vjust = 0.5),
#                 strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
#           labs(x = "", y = "",
#               fill = "-log10[PValue.adjust]") +
#           ggsave(filename = paste(fig.path, "/DEGs/geneset_aracyc.png", sep = ""), 
#            dpi = 600, 
#            device = "png", 
#            height = 21, 
#            width = 14, limitsize = F, bg = "transparent")

# Target <- "MYC"
# Set the target marker genes and check the LOGFC profile
target <- c("AT1G32640", "AT5G46760", "AT4G17880", "AT1G18710", "AT1G74430", "AT5G24780")
names(target) <- c("MYC2", "MYC3", "MYC4", "MYB47", "MYB95", "VSP1")

df_mycs <- df %>% filter(gene_ids %in% target) %>%
mutate(gene_ids = factor(gene_ids, levels = target, label = names(target)))

# MYC profile
df_mycs %>% 
            ggplot(aes(y = Genotype, x = logFC)) +
            ggtitle("Targetted Expression Profile") +
            geom_vline(xintercept = 0.0, lty = "solid", colour = "darkgrey", lwd = 2) +
            geom_point(alpha = 0.8,
              stat = "identity", 
              # position = position_dodge(width = 0.85), 
              na.rm = FALSE
              ) +
            facet_grid(gene_ids ~ .,
               switch = "both", space = "free", scale = "free",
               labeller = labeller(label.idx)) +
            # geom_boxplot(alpha = 0.6, outlier.alpha = 0, fill = NA, aes(colour = treatment), na.rm = FALSE) +
            scale_colour_manual(values = genotype$col.idx, 
              label = genotype$name) +
            theme_AKB +
            theme(
                  axis.text.y = element_text(hjust = 0.8, vjust = 0.5),
                  strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                  legend.text = element_text(size = 6)) +
            labs(x = "", y = "",
                colour = "logFoldChange") +
            ggsave(paste(fig.path, "/genewise/boxplot_(MYCs_MYBs)_target.png", sep = ""), 
                   dpi = 600, 
                   device = "png", 
                   height = 10, 
                   width = 7, limitsize = F, bg = "transparent")

# Targetted analysis and statistics
load(paste(out.path, "/log.normalised.matrix.Rdata", sep = ""))
genelist <- read.delim2(file = paste0("./data/genelist_target.txt"), sep = ",", as.is = TRUE, header = TRUE)

# Make data frame
df.target <- df %>% 
    mutate(gene_ids = as.factor(as.character(gene_ids)),
      symbol = as.factor(paste(gene_ids, as.character(logFC_P.annotated$symbol[match(.$gene_ids, row.names(logFC_P.annotated))]), sep = "::")),
      counts = ifelse(!is.finite(log2(counts)), NA, log2(counts))
      ) %>%
    # filter(gene_ids %in% target) %>%
    data.frame(.)


# logFC_P.annotated$symbol[str_detect(logFC_P.annotated$symbol, "TSA1")]
# genelist$gene_ids <- df.target$gene_ids[match(genelist$gene_symbol, as.character(df.target$symbol))]
genelist$gene_ids <- mclapply(genelist$gene_symbol, function(x) {
  idx <- str_detect(df.target$symbol, x)
  if(sum(idx) > 0) {
    ret <- paste(unique(df.target$gene_ids[idx]), collapse = ";")
    } else ret <- NA
  return(ret)
  }, mc.cores = 4) %>% 
unlist(.)

idx <- match(as.character(df.target$gene_ids), as.character(genelist$gene_ids))
df.target <- cbind.data.frame(df.target, 
  genelist[idx, c("gene_symbol", "GLSpathway", "rank")])

# Mark TSA1 and bglu21
df.target$gene_symbol[str_detect(df.target$symbol, "TSA1")] <- genelist$gene_symbol[str_detect(genelist$gene_symbol, "TSA1")]
df.target$GLSpathway[str_detect(df.target$symbol, "TSA1")] <- genelist$GLSpathway[str_detect(genelist$gene_symbol, "TSA1")]
df.target$rank[str_detect(df.target$symbol, "TSA1")] <- genelist$rank[str_detect(genelist$gene_symbol, "TSA1")]


idx <- match(row.names(logFC_P.annotated), as.character(genelist$gene_ids))
lfc_tabs <- cbind.data.frame(logFC_P.annotated, 
  genelist[idx, c("gene_symbol", "GLSpathway", "rank")])

lfc_tabs$gene_symbol[str_detect(lfc_tabs$symbol, "TSA1")] <- genelist$gene_symbol[str_detect(genelist$gene_symbol, "TSA1")]
lfc_tabs$GLSpathway[str_detect(lfc_tabs$symbol, "TSA1")] <- genelist$GLSpathway[str_detect(genelist$gene_symbol, "TSA1")]
lfc_tabs$rank[str_detect(lfc_tabs$symbol, "TSA1")] <- genelist$rank[str_detect(genelist$gene_symbol, "TSA1")]


idx <- str_detect(colnames(lfc_tabs), "^genotype")
df_lfc <- lfc_tabs[,!idx] %>% na.omit(.) %>%
add_column(id = row.names(.), .before = 1) %>%
gather(key = "interaction", value = "vals", convert = FALSE, -id, -gene_name, -go_id, -enzyme, -symbol, -kegg_path, -aracyc, -gene_symbol, -GLSpathway, -rank) %>%
separate(interaction, into = c("x", "Genotype", "x1", "wide"), sep = "_", convert = FALSE) %>%
spread(key = wide, value = vals, convert = FALSE, fill = NA) %>%
dplyr::select(-x,-x1) %>%
mutate(Genotype = factor(Genotype, levels= genotype$short),
  sig = ifelse(PValue <= alpha, sign(logFC) * 1, 0),
  logFC = logFC) %>%
data.frame(.)

# Genesets for GLS pathway, ERB pathway and other transcription associations

# df.target <- df %>% 
#     mutate(gene_ids = as.factor(as.character(gene_ids)),
#       symbol = as.factor(paste(gene_ids, as.character(logFC_P.annotated$symbol[match(.$gene_ids, row.names(logFC_P.annotated))]), sep = "::")),
#       gene_name = as.factor(paste(gene_ids, as.character(logFC_P.annotated$gene_name[match(.$gene_ids, row.names(logFC_P.annotated))]), sep = "::")),
#       counts = log2(1+counts)
#       ) %>%
#    mutate(glucosinolate_genes = str_detect(gene_name, "glucosinolate"),
#     myb_transcription = str_detect(gene_name, "transcription factor")) %>%
#    filter(glucosinolate_genes == TRUE) %>%
#     data.frame(.)

# Conduct gene by genewise statistics and plot
df_target <- df.target %>% na.omit(.) 
fit_list <- mclapply(c(as.character(unique(df_target$gene_ids))), function(x){

        # x <- "AT1G16400"
        temp <- df.target %>% filter(gene_ids == x)
        sym <- paste(paste0(unique(as.character(temp$rank)), "__", unique(as.character(temp$gene_symbol))), unique(as.character(temp$GLSpathway)), sep = ";")
        temp %>% 
            ggplot(aes(x = genotype, 
                          y = ifelse(!is.finite(log2(counts)), NA, log2(counts)))) +
            ggtitle(paste0(sym)) +
            geom_jitter(alpha = 0.8, aes(colour = treatment), 
              stat = "identity", 
              position = position_dodge(width = 0.85), na.rm = FALSE
              ) +
            geom_boxplot(alpha = 0.6, outlier.alpha = 0, fill = NA, aes(colour = treatment), na.rm = FALSE) +
            scale_colour_manual(values = treatment$col.idx, 
              label = treatment$name) +
            theme_AKB +
            theme(
                  axis.text.y = element_text(hjust = 0.8, vjust = 0.5),
                  strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                  legend.text = element_text(size = 6)) +
            labs(x = "", y = "",
                colour = "log2[counts]") +
            ggsave(paste(fig.path, "/genewise/boxplot_(",x,")_target.png", sep = ""), 
                   dpi = 600, 
                   device = "png", 
                   height = 10, 
                   width = 7, limitsize = F, bg = "transparent")

        mod <- aov(lm(counts ~ 0 + group + random, temp))
        fit <-  broom::tidy(TukeyHSD(mod)) %>% 
                filter(term == "group") %>%
                separate(comparison, into = c("A", "B"), convert = F, sep = "-") %>%
                separate(A, into = c("genotypeA", "treatmentA"), convert = F, sep = "_") %>%
                separate(B, into = c("genotypeB", "treatmentB"), convert = F, sep = "_") %>%
                mutate(significance = ifelse(adj.p.value < 0.01, "**", ifelse(adj.p.value < 0.05, "*", "NS")),
                  term = x, gene_symbol = sym) %>%
                arrange(significance) %>%
                data.frame(.)

        # temp %>% 
        # mutate(mx = counts + sd(counts),
        #  min = counts - sd(counts)) %>%
        #     ggplot(aes(x = genotype, 
        #                   y = counts)) +
        #     ggtitle(paste0(sym)) +
        #     geom_bar(alpha = 1, stat = "identity", position = position_dodge(),
        #      aes(fill = treatment)) +
        #     geom_errorbar(aes(x = genotype, ymin = min, ymax = mx, colour = treatment), alpha = 0.8, width = 0.4) +
        #     scale_fill_manual(values = treatment$col.idx, 
        #       label = treatment$name) +
        #     cale_colour_manual(values = treatment$col.idx, 
        #       label = treatment$name) +
        #     theme_AKB +
        #     theme(
        #           axis.text.y = element_text(hjust = 0.8, vjust = 0.5),
        #           strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        #           legend.text = element_text(size = 6)) +
        #     labs(x = "", y = "",
        #         colour = "log2[counts]") +
        #     ggsave(paste(fig.path, "/genewise/barplot_(",x,")_target.png", sep = ""), 
        #            dpi = 600, 
        #            device = "png", 
        #            height = 10, 
        #            width = 7, limitsize = F, bg = "transparent")


        return(fit)

}, mc.cores = 12)


# Bind it into a dataframe
res <- do.call(rbind.data.frame, fit_list) %>% arrange(significance)
res_treatment <- res %>% 
filter(treatmentA != treatmentB, genotypeA == genotypeB) %>%
mutate(Genotype = factor(genotypeA, levels = genotype$short),
  sig = ifelse(adj.p.value <= 0.05, sign(estimate)*1, 0)
  ) %>% 
data.frame(.)


levels <- lapply(unique(genelist$gene_symbol), function(x)
  paste(min(genelist$rank):max(genelist$rank), x, sep = "__")) %>% 
  unlist(., use.names = FALSE)

# Make hmap for the expression
df.target %>% 
na.omit(.) %>%
mutate(gene_symbol = factor(paste(rank,gene_symbol, sep = "__"), levels = levels), 
  GLSpathway = factor(GLSpathway, levels = c("ERB", "aliphatic", "indolic", "both", "JA response"))) %>%
arrange(rank) %>%
    ggplot(aes(x = sample_name, 
                  y = gene_symbol, 
                  fill = log2(1+counts))) +
    geom_raster(alpha = 1) +
    facet_grid(GLSpathway ~ genotype + treatment,
               switch = "both", space = "free", scale = "free",
               labeller = labeller(label.idx)) +
    scale_fill_gradient(na.value = "darkgrey", 
                        low = "black", 
                       high = "yellow") +
    theme_AKB +
    theme(panel.spacing = unit(0.2, "lines"),
      axis.text.x = element_blank(),
          axis.text.y = element_text(hjust = 0.8, vjust = 0.5),
          strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 6)) +
    labs(x = "", y = "",
        fill = "log2[counts]") +
    ggsave(paste(fig.path, "/EDA_heatmap_treatment_target.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           height = 10, 
           width = 7, limitsize = F, bg = "transparent")

# Make hmap for the expression
df_lfc %>% 
mutate(gene_symbol = factor(paste(rank, gene_symbol, sep = "__"), levels = levels), 
  GLSpathway = factor(GLSpathway, levels = c("ERB", "aliphatic", "indolic", "both", "JA response")),
  sig = factor(sig, levels = c(-1,1,0))
  ) %>%
arrange(desc(rank)) %>%
    ggplot(aes(x = Genotype, 
                  y = gene_symbol, 
                  fill = logFC)
    ) +
    geom_raster(alpha = 1) +
    geom_tile(alpha = 0.6, fill = NA, aes(colour = sig), 
      size = 2, width = 0.95, height = 0.95) +
    facet_grid(GLSpathway ~ .,
               switch = "both", space = "free", scale = "free",
               labeller = labeller(label.idx)) +
    scale_fill_gradient2(na.value = "darkgrey", 
      mid = "white", 
      midpoint = 0.0,
      low = "darkmagenta", 
      high = "darkgreen") +
    scale_colour_manual(values = c(`1` = "red", `-1` = "blue", `0` = "transparent"), 
      guide = FALSE, na.value = NA, 
      label = c(`1` = "Up in JA", `-1` = "Up in DW", `0` = "transparent")
      ) +
    theme_AKB +
    theme(panel.spacing = unit(0.2, "lines"),
      axis.text.x = element_blank(),
          axis.text.y = element_text(hjust = 0.8, vjust = 0.5),
          strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 6)) +
    labs(x = "", y = "",
        fill = "log2[foldchange]") +
    ggsave(paste(fig.path, "/LFC_heatmap_treatment_target.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           height = 10, 
           width = 7, limitsize = F, bg = "transparent")

res_treatment %>% 
separate(gene_symbol, into = c("gene_symbol", "GLSpathway"), sep = ";", convert = FALSE) %>%
mutate(gene_symbol = factor(gene_symbol, levels = levels), 
  GLSpathway = factor(GLSpathway, levels = c("ERB", "aliphatic", "indolic", "both", "JA response")),
  sig = factor(sig, levels = c(-1,1,0)), 
  estimate = log2(estimate+1)
  ) %>%
    ggplot(aes(x = Genotype, 
                  y = gene_symbol, 
                  fill = estimate)
    ) +
    geom_raster(alpha = 1) +
    geom_tile(alpha = 0.6, fill = NA, aes(colour = sig), 
      size = 2, width = 0.95, height = 0.95) +
    facet_grid(GLSpathway ~ .,
               switch = "both", space = "free", scale = "free",
               labeller = labeller(label.idx)) +
    scale_fill_gradient2(na.value = "darkgrey", 
      mid = "white", 
      midpoint = 0.0,
      low = "darkmagenta", 
      high = "darkgreen") +
    scale_colour_manual(values = c(`1` = "red", `-1` = "blue", `0` = "transparent"), 
      guide = FALSE, na.value = NA, 
      label = c(`1` = "Up in JA", `-1` = "Up in DW", `0` = "transparent")
      ) +
    theme_AKB +
    theme(panel.spacing = unit(0.2, "lines"),
      axis.text.x = element_blank(),
          axis.text.y = element_text(hjust = 0.8, vjust = 0.5),
          strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 6)) +
    labs(x = "", y = "",
        fill = "log2(Mean Estimate)") +
    ggsave(paste(fig.path, "/Mean_heatmap_treatment_target.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           height = 10, 
           width = 7, limitsize = F, bg = "transparent")



# Annotate
# fit <- do.call(rbind.data.frame, fit_list) %>% arrange(significance)
# fit <- cbind.data.frame(fit, logFC_P.annotated[match(fit$term, row.names(logFC_P.annotated)), c("symbol", "kegg_path", "aracyc")])
# row.names(fit) <- NULL

write.table(res, 
            paste(out.path, "/log.normalised.statistics.thsd_pHOC.targetted.txt", 
                  sep = ""), sep = "\t", quote = F)
      
message("DONE")

# Write output
save(list = c("logFC_P.annotated"),
     file = paste(out.path, "/gset_analysis/objects/deseq-lfc_annotated.Rdata", sep = ""))
write.table(logFC_P.annotated, 
            paste(out.path, "/deseq-lfc-p_annotated.txt", 
                  sep = ""), sep = "\t", quote = F)

# END OF SCRIPT
sessionInfo()