#!/bin/bash

# Use a logger using date function

# Set path
path_data="/klaster/scratch/CAMDA/CAMDA2019/MetagenomicForensics/REPOSITORY/RNAseq/";
data_output="/klaster/work/abasak/git_repo_myb/rna_seq/data/";
sample_name=$(cat '/home/abasak/git_repo_myb/rna_seq/data/samplenames.txt');


# Prepare Reference
ref_list=('cdna' 'cds' 'ncrna')

echo "Preparing reference ..."
for r in "${ref_list[@]}"
do
	echo "Reference set for $r ...";
	
	zcat -rf $path_data/TAIR10db/ref_seq/TAIR.10/$r/Arabidopsis_thaliana.TAIR10.$r.all.fa.gz \
	&>> $path_data/temp/Arabidopsis_thaliana.TAIR10.$r.all.fa;

	rsem-prepare-reference --gtf $path_data/TAIR10db/ref_seq/TAIR.10/gtf/Arabidopsis_thaliana.TAIR10.42.gtf \
	--bowtie2 \
	$path_data/temp/Arabidopsis_thaliana.TAIR10.$r.all.fa \
	$path_data/TAIR10db/RSEM_ref/arath_$r;

	rm -rf $path_data/temp/*;

done;
echo "Reference Build Complete !!"


echo "Reference Build SALMON Pipeline !!"
salmon index -t $path_data/TAIR10db/ref_seq/TAIR.10/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz -i athal_index