# Script to make read counts table from SALMON output
# Clean up
rm(list = ls())

options(warn = 1,
       mc.cores = 8)

# Load necessary packages
pkgs <- c("tidyverse", "Rsamtools", "tximport", "GenomicAlignments", "GenomicFeatures", "BiocParallel", "TxDb.Athaliana.BioMart.plantsmart28")

lapply(pkgs, require, character.only = TRUE)

# Set path
path <- "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
metadata <- read.table(paste(data.files, "metadata.txt", sep = ""),
                       header = TRUE, as.is = TRUE)
#metadata$sample_name <- str_replace(metadata$sample_name, "_WT_JA_4_28", "WT_JA_4_28") # String parsing

squant.files <- list.files(salmon.path,
                        full.names = T, 
                        include.dirs = T, 
                        ignore.case = T, 
                        recursive = T, 
                        pattern = "quant.sf")

idx <- sapply(unique(metadata$sample_name), function(x) grep(squant.files, 
                                                             pattern = paste(x, "salmon", sep = "_"), 
                                                             value = TRUE), 
              simplify = TRUE, USE.NAMES = T) %>% unlist(.)

# Run Genomic alignments and build a TxDB object from the database directly
k <- keys(TxDb.Athaliana.BioMart.plantsmart28, keytype = "TXNAME")
txt <- select(TxDb.Athaliana.BioMart.plantsmart28, k, "GENEID", "TXNAME")
ex <- exonsBy(TxDb.Athaliana.BioMart.plantsmart28, "gene")
# idx <- idx[!str_detect(idx, filter.out)]

tximport_obj <- tximport(idx, type = "salmon", tx2gene = txt, txOut = TRUE, 
                         varReduce = TRUE)
tximport_sum <- summarizeToGene(tximport_obj, txt)

flag <- all.equal(colSums(as.data.frame(tximport_obj$counts)), 
                  colSums(as.data.frame(tximport_sum$counts)))

if(flag == TRUE){ 
    message("Transcript build successful..!!")
    }else message("ERROR Please check..!!")

# Write the output as table
write.table(as.data.frame(tximport_obj$counts), 
            paste(out.path, "/summarised_exp/edata_rs.txt", sep = ""),
            quote = FALSE)

write.table(as.data.frame(tximport_sum$counts), 
            paste(out.path, "/summarised_exp/edata_rs_genelevel.txt", sep = ""),
            quote = FALSE)

#write.table(as.data.frame(rowRanges(ga_obj)$annotation), paste(out.path, "/summarised_exp/fdata_rs.txt", sep = ""))
save(list = c("tximport_obj"), 
     file = paste(out.path, "/summarised_exp/salmon_summarized_experiments.Rdata", sep = ""))

save(list = c("tximport_sum"), 
     file = paste(out.path, "/summarised_exp/salmon_summarized_experiments_genelevel.Rdata", sep = ""))

# END OF SCRIPT
sessionInfo()
