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
idy <- str_detect(colnames(logFC_P.annotated), "^treatment") & str_detect(colnames(logFC_P.annotated), "_logFC$")
df <- cbind.data.frame(logFC_P.annotated[,idy], symbol = logFC_P.annotated$symbol)

# Mine MYBs and MYCs
idx <- str_detect(df$symbol, "myb|MYB|myc|MYC")
df <- df[idx,]

d <- 1 - cor(t(df[, !colnames(df) == "symbol"]))
hc <- hclust(as.dist(d), "ward.D")
myc_myb_lfc_clusters <- row.names(d)[hc$order]

cor_df <- as.data.frame(d) %>%
add_column(source = row.names(.), .before = 1) %>%
gather(key = "sink", value = "pcc", convert = FALSE, -source) %>%
mutate(myb = df[.$source,"symbol"],
  myc = df[.$sink,"symbol"],
  pcc = 1 - pcc, 
  source = factor(source, levels = myc_myb_lfc_clusters), 
  sink = factor(sink, levels = myc_myb_lfc_clusters)
  )

# cor_df$trace <- NA
# cor_df$trace[str_detect(cor_df$myb, "MYB|myb")] <- "MYBs" 
# cor_df$trace[str_detect(cor_df$myc, "MYC|myc")] <- "MYCs" 
idx <- str_detect(cor_df$myb, "MYB|myb") & str_detect(cor_df$myc, "MYC|myc")
cor_df[idx,] %>%
mutate(trace = str_detect(myb, "MYB47|myb47|MYB95|myb95") & str_detect(myc, "MYC2|myc2|MYC3|myc3|MYC4|myc4")) %>%
ggplot(aes(y = source, 
                  x = myc)
    ) +
    geom_raster(alpha = 1, aes(fill = pcc)) +
    geom_tile(alpha = 0.6, fill = NA, aes(colour = trace), 
      size = 0.8, width = 0.95, height = 0.95) +
    # facet_grid(. ~ trace,
    #            switch = "both", 
    #            space = "free", 
    #            scale = "free") +
    scale_fill_gradient2(na.value = "darkgrey", 
      mid = "white", 
      midpoint = 0.0,
      low = "darkmagenta", 
      high = "darkgreen", limits = c(-1,1), 
      breaks = c(-1,0,1), 
      label = c("-1", "0", "1")) +
    scale_colour_manual(values = c("transparent", "black"), 
      guide = FALSE, na.value = NA
      ) +
    theme_AKB +
    theme(panel.spacing = unit(0.2, "lines"),
      axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 0.9, angle = 90),
          axis.text.y = element_text(hjust = 0.9, vjust = 0.5, size = 8),
          # strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          # strip.text.y = element_text(angle = 180, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 8)) +
    labs(x = "", y = "",
        fill = "PCC") +
    ggsave(paste(fig.path, "/myb_myc_coexp.png", sep = ""), 
           dpi = 600, 
           device = "png", 
           height = 14, 
           width = 4, limitsize = F, bg = "transparent")

# END OF SCRIPT
sessionInfo()