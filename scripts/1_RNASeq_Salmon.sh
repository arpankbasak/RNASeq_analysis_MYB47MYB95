#!/bin/bash

# Use a logger using date function

# Set path
path_data="/klaster/scratch/CAMDA/CAMDA2019/MetagenomicForensics/REPOSITORY/RNAseq";
data_output="/klaster/work/abasak/git_repo_myb/rna_seq/data";
sample_name=$(cat '/klaster/work/abasak/git_repo_myb/rna_seq/data/samplenames.txt');

# Cleanup
rm -rf $path_data/aligned_salmon/*

# ---
date ;
echo "Initializing RSEM for RNA Expression profiling ...";

for f in $sample_name
do
#     echo "Initiating QC for sequences recovered  ${f} ...";
#     fastqc -o $path_data/qc_before_trimmming/ --noextract \
#     $path_data/raw_data/${f}_1.fastq.gz $path_data/raw_data/${f}_2.fastq.gz ;
      
    echo "Initiating TRIM GALORE for removing adapters ${f} ..."
    trim_galore --paired $path_data/raw_data/${f}_1.fastq.gz $path_data/raw_data/${f}_2.fastq.gz \
    --core 8 \
    --gzip \
    --illumina \
    --phred33 \
    -o $path_data/trimmed_reads ;
    
    echo "Initiating SALMON Pipeline for sample ${f} ..."
    salmon quant -i $path_data/TAIR10db/ref_seq/TAIR.10/athal_index \
    -l ISR \
    -1 $path_data/trimmed_reads/${f}_1_val_1.fq.gz \
    -2 $path_data/trimmed_reads/${f}_2_val_2.fq.gz \
    -p 8 \
    --validateMappings \
    -o $path_data/aligned_salmon/${f}_salmon ;

done;

echo "Successfully Completed SALMON pipeline";
date ;

