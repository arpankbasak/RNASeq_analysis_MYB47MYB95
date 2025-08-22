rm(list = ls())

options(warn = 1,
       mc.cores = 6)

pkgs <- c("tidyverse", "fgsea", "GO.db", "parallel")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
load(paste(out.path, "/DEG-list_deseq-treatment.Rdata", sep = ""))
message("Saving gene lists ...")
sink(paste(paste(out.path, "/gset_analysis/deseq-upregulated_treatment_gene_list.txt", sep = "")))
upregulated
sink(paste(paste(out.path, "/gset_analysis/deseq-downregulated_treatment_gene_list.txt", sep = "")))
downregulated
sink()
sink()

message("Building annotation object using Biomart for upregulated genes ...")
# Annotating Upregulated genes from biomart
upregulated <- unique(unlist(upregulated, use.names = F))
bm_res_upregulated <- mclapply(upregulated, function(x){
    
        temp <- call_biomart(x)
        return(temp)
})
#names(bm_res_upregulated) <- names(upregulated)
save(list = "bm_res_upregulated", 
     file = paste(out.path, "gset_analysis/objects/deseq-biomart_treatment_up.Rdata", sep = ""))
rm(list = "upregulated")

message("Building annotation object using Biomart for downregulated genes ...")
# Annotating downregulated genes from Biomart
downregulated <- unique(unlist(downregulated, use.names = F))
bm_res_downregulated <- mclapply(downregulated, function(x){
    
        temp <- call_biomart(x)
        return(temp)
})
#names(bm_res_downregulated) <- names(downregulated)
save(list = "bm_res_downregulated", 
     file = paste(out.path, "gset_analysis/objects/deseq-biomart_treatment_down.Rdata", sep = ""))
rm(list = "downregulated")
message("DONE")

# END OF SCRIPT
sessionInfo()