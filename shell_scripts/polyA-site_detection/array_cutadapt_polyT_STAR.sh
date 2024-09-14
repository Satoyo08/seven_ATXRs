#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd

# array_cutadapt_polyT_STAR.sh
module load singularity
alias cutadapt='singularity exec /usr/local/biotools/c/cutadapt:1.15--py36_0 cutadapt'
alias STAR="singularity exec /usr/local/biotools/s/star:2.5.2b--0 STAR"

# file handle
# usage: ls *fastq.gz> fastq_list_R1 
array=(`cat fastq_list_R1`)
in_pair1=${array[$SGE_TASK_ID]}

# cutadapt , polyT removal parameters
mismatch=0.14
min_matched_length=6
min_length_output=20

outfile="trimmed_polyT.e"$mismatch"_minlen"$min_length_output"_"$in_pair1

cutadapt \
	-g "XT{51}" \
	--trimmed-only -e $mismatch -O $min_matched_length --minimum-length $min_length_output \
	-o $outfile \
	$in_pair1
	
# somehow, for single end, mapping doesnt work without prior decompression...
gunzip $outfile 

STAR --outSAMtype BAM SortedByCoordinate \
     --genomeDir   /home/satoyo08/refs/Arabidopsis_STAR_araport \
     --readFilesIn  $(echo $outfile | sed -e 's/.gz//')  \
     --genomeLoad NoSharedMemory \
     --outFilterMultimapNmax 3 \
     --outFileNamePrefix "trimmed_polyT"$(echo $in_pair1 | sed -e 's/.fastq.gz//') 