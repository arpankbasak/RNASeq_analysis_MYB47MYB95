# Script to parse in the KEGG map for differential genes compared to WT and control
# @Arpan Kumar Basak

# Cleanup
rm(list = ls())

options(warn = 1,
       mc.cores = 24)

pkgs <- c("tidyverse", "pathview", "org.At.tair.db", "parallel", "KEGGREST")
lapply(pkgs, require, character.only = TRUE)

# Set path
path = "/klaster/work/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)

# Raad in data
load(paste(out.path, "gset_analysis/objects/edeger-GO-Pathway_treatment.Rdata", sep = ""))
load(file = paste(out.path, "edger_objects.plotting.Rdata", sep = ""))

# Feed in the LFC table and LFC as input for pathview -- incase of thresholding out the non significant genes
idx <- str_detect(colnames(logFC_P), "^treatment") & str_detect(colnames(logFC_P), "logFC$")
logFC_P <- logFC_P[, idx]

# Make empty list or placeholder
pv.out <- list()

# Change working directory download the dedicated maps
setwd(paste(path,fig.path,"/KEGG_MAPs/treatment/", sep = ""))

# Build a pathwview object to run the pathview pipeline
pathview_obj <- mclapply(1:length(path_df$PATH), function(m){
    
    #m =3;i = 1
    x <- path_df$PATH[m]
    k <- path_df$kegg_path_id[m]
    j <- c(KEGGREST::keggGet(as.character(k))[[1]]$NAME[1], use.names = FALSE)
    message("Downloading map for [", j, "] pathway ...", sep = "")

    
    for(i in 1:length(colnames(logFC_P))){
        
        temp <- c(logFC_P[, i], use.names = TRUE)
        #names(temp) <- str_replace(row.names(logFC_P), "T", "t")
        names(temp) <- row.names(logFC_P)
        pv.out[[i]] <- pathview(gene.data = temp, 
                           pathway.id = x, 
                           species = "ath", 
                                gene.annotpkg="org.At.tair.db", 
                                gene.idtype = "TAIR",
                           kegg.native = TRUE,
                           out.suffix = paste(colnames(logFC_P)[i], 
                                              "_", str_replace_all(j, "\\s+", "_"), sep = ""), 
                           kegg.dir = paste(getwd(), "/archive", sep = ""),
                           low = gradient.low, 
                           mid = gradient.mid, 
                           high = gradient.high, 
                           na.col = gradient.na)
        
    }
   names(pv.out) <- colnames(logFC_P)
    
    return(pv.out)
    message("DONE")
}, mc.cores = 24)

message("DONE")

# END OF SCRIPT
sessionInfo()