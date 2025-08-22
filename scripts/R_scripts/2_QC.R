rm(list = ls())

options(warn = 1,
       mc.cores = 8)

pkgs <- c("tidyverse", "rsamtools", "Rsubread", "biocparallel")
mclapply(pkgs, require, characters.only = TRUE)

# Set path
path = "/home/abasak/git_repo_myb/rna_seq/"
setwd(path)
param <- "./scripts/R_scripts/parameters.R"
source(param)



# Raad in data
bam.files <- list.files(paste(data.files, "bam_files/", sep = ""),
                        full.names = T, 
                        include.dirs = T, 
                        ignore.case = T, 
                        recursive = T, 
                        pattern = ".bam"))

# testing
bam.files <- bam.files[1]

metadata <- read.table(paste(data.files, "metadata.txt", sep = ""),
                       header = TRUE, as.is = TRUE)

# Run RSUBREAD
rs_obj <- Rsubread::featureCounts(bam.files,
                        annot.ext=gff.path,
                        isGTFAnnotationFile=TRUE,
                        isPaired=TRUE)

# Write the output as table
write.table(as.data.frame(rs_obj$counts), paste(out.files, "./rsubread/edata_rs.txt", sep = ""))
write.table(as.data.frame(rs_obj$annotation), paste(out.files, "./rsubread/fdata_rs.txt", sep = ""))

# Match conalmes and rename columns with sample data
colnames(rs_obj$counts) <- str_replace(colnames(rs_obj$counts), pattern = ".bam", replacement = "")
if(sum(match(colnames(rs_obj$counts), unique(metadata$treatment))) != 0) {
    message("Colnames unique")
    
    # Save as objects
    save(list = "rs_obj", paste(out.path, "rsubread/rs_featurecounts.Rdata", sep = ""))

}else message("ERROR!!:: Colnames not unique")

# Save as objects
save(list = "rs_obj", paste(out.path, "rsubread/rs_featurecounts.Rdata", sep = ""))

