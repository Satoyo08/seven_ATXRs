#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd

# array_cutadapt_rm_adapter.sh
module load singularity
alias cutadapt='singularity exec /usr/local/biotools/c/cutadapt:1.15--py36_0 cutadapt'

# file handle
# usage: ls *1.fq.gz > fastq_list_R1 
array=(`cat fastq_list_R1`)
in_pair1=${array[$SGE_TASK_ID]}
in_pair2=$(echo $in_pair1 | sed -e 's/_1\./_2\./')


# cutadapt, trim the adapter
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o "trimmed_"$in_pair1 -p "trimmed_"$in_pair2 \
    --minimum-length 4 \
    $in_pair1 $in_pair2