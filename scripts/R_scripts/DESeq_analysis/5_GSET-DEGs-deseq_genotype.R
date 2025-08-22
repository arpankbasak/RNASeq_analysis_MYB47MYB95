rm(list = ls())

options(warn = 1, mc.cores = 24)

pkgs <- c("tidyverse", "fgsea", "GO.db", "parallel")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
load(paste(out.path, "/DEG-list_deseq-genotype.Rdata", sep = ""))
message("Saving gene lists ...")
sink(paste(paste(out.path, "/gset_analysis/deseq-upregulated_genotype_gene_list.txt", sep = "")))
upregulated
sink(paste(paste(out.path, "/gset_analysis/deseq-downregulated_genotype_gene_list.txt", sep = "")))
downregulated
sink()

message("Building annotation object using Biomart for upregulated genes ...")
upregulated <- unique(unlist(upregulated, use.names = F))
# Annotating Upregulated genes from biomart
bm_res_upregulated <- mclapply(as.list(upregulated), function(x){
    
        temp <- call_biomart(x)
        return(temp)
})

#names(bm_res_upregulated) <- names(upregulated)
save(list = "bm_res_upregulated", 
     file = paste(out.path, "gset_analysis/objects/edeger-biomart_genotype_up.Rdata", sep = ""))
rm(list = c("upregulated", "bm_res_upregulated"))

message("Building annotation object using Biomart for downregulated genes ...")

# Annotating downregulated genes from Biomart
downregulated <- unique(unlist(downregulated, use.names = F))
bm_res_downregulated <- lapply(as.list(downregulated), function(x){
    
        temp <- call_biomart(x)
        return(temp)
})
#names(bm_res_downregulated) <- names(downregulated)
save(list = "bm_res_downregulated", 
     file = paste(out.path, "gset_analysis/objects/edeger-biomart_genotype_down.Rdata", sep = ""))
rm(list = c("downregulated", "bm_res_downregulated"))
message("DONE")

# END OF SCRIPT
sessionInfo()