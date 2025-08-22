#!/bin/bash

# Use a logger using date function

# Set path
path_data="/klaster/scratch/CAMDA/CAMDA2019/MetagenomicForensics/REPOSITORY/RNAseq";
data_output="/home/abasak/git_repo_myb/rna_seq/data";
sample_name=$(cat '/home/abasak/git_repo_myb/rna_seq/data/samplenames.txt');
script_path="/klaster/work/abasak/git_repo_myb/rna_seq/scripts/R_scripts/"

# ---
echo "Initializing DESEq2 protocol ...";

Rscript $script_path/1_EDA.R ;
Rscript $script_path/DESeq_analysis/4_DEG-deseq-genotype_plot_canvas.R ;
Rscript $script_path/DESeq_analysis/4_DEG-deseq-treatment_plot_canvas.R ;
Rscript $script_path/DESeq_analysis/5_Annotation-DEGs-deseq_global.R;

rm /klaster/work/abasak/git_repo_myb/rna_seq/figures/DEGs/*.log ;

echo "Successfully Completed DESeq2 analysis pipeline";