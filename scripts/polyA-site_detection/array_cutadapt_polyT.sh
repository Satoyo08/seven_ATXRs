#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd

# array_cutadapt_polyT.sh

module load singularity
alias cutadapt='singularity exec /usr/local/biotools/c/cutadapt:1.15--py36_0 cutadapt'

# file handle
# usage: ls *1.fq.gz > fastq_list_R1 
array=(`cat fastq_list_R1`)
in_pair1=${array[$SGE_TASK_ID]}
in_pair2=$(echo $in_pair1 | sed -e 's/_1\./_2\./')

# cutadapt , polyT removal parameters
mismatch=0.14
min_matched_length=6
min_length_output=4

cutadapt \
	-g "XT{151}" \
	---trimmed-only -e $mismatch -O $min_matched_length --minimum-length $min_length_output \
	-o "trimmed_polyT.e"$mismatch"_minlen"$min_length_output"_"$in_pair1 -p "trimmed_polyT.e"$mismatch"_minlen"$min_length_output"_"$in_pair2 \
	"trimmed_"$in_pair1 "trimmed_"$in_pair2