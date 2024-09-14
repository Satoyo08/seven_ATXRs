#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd

# star_polyT.sh
module load singularity
alias STAR="singularity exec /usr/local/biotools/s/star:2.5.2b--0 STAR"


# file handle
# usage: ls *1.fq.gz > fastq_list_R1 
array=(`cat fastq_list_R1`)
in_pair1=${array[$SGE_TASK_ID]}
in_pair2=$(echo $in_pair1 | sed -e 's/_1\./_2\./')

mismatch=0.14
min_matched_length=6
min_length_output=4

trimmedfqgz_1="trimmed_polyT.e"$mismatch"_minlen"$min_length_output"_"$in_pair1
trimmedfqgz_2="trimmed_polyT.e"$mismatch"_minlen"$min_length_output"_"$in_pair2
gunzip $trimmedfqgz_1 $trimmedfqgz_2

STAR --outSAMtype BAM SortedByCoordinate \
     --genomeDir   /home/satoyo08/refs/Arabidopsis_STAR_araport \
     --readFilesIn $(echo $trimmedfqgz_1 | sed -e 's/.gz//') $(echo $trimmedfqgz_2 | sed -e 's/.gz//') \
     --genomeLoad NoSharedMemory \
     --outFilterMultimapNmax 3 \
     --outFileNamePrefix "trimmed_polyT_"$(echo $in_pair1 | sed -e 's/_1.fq.gz//') 