#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd
#$ -pe def_slot 4

# prep ; ls *_1.fastq.gz > fastq_list_R1
# usage qsub -l short -t 1-16:1 
export MALLOC_ARENA_MAX=2
module load singularity
adaptor=~/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa
alias trim='singularity exec /usr/local/biotools/t/trimmomatic:0.36--5 trimmomatic -Xmx2g'
array=(`cat fastq_list_R1`)
pair1=${array[$SGE_TASK_ID]}
pair2=$(echo $pair1 | sed -e 's/_1\./_2\./')
trim PE -threads 4 -phred33 -trimlog $pair1"trimming_log.txt" $pair1 $pair2 'paired_trimmed_'$pair1 'unpaired_trimmed_'$pair1 'paired_trimmed_'$pair2 'unpaired_trimmed_'$pair2 ILLUMINACLIP:$adaptor:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36
