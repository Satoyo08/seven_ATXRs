#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd

#usage qsub -t 1-12 -tc 3 -l s_vmem=32G,mem_req=32G ~/shells/array_jobs/general_fastq_to_bed_array.sh
module load singularity
array=(`cat fastq_list_R1`)
pair1=${array[$SGE_TASK_ID]}
pair2=$(echo $pair1 | sed -e 's/_1\./_2\./')
a1=$(echo $pair1 | sed -e 's/1\.fastq.gz//')
gunzip 'paired_trimmed_'$a1'1.fastq.gz' 'paired_trimmed_'$a1'2.fastq.gz'
printf "this_is\t$a1\n"

singularity exec /usr/local/biotools/b/bowtie2:2.3.3.1--py36pl5.22.0_0 bowtie2 -x ~/ref_tair10_bowtie2 -1 'paired_trimmed_'$a1'1.fastq' -2 'paired_trimmed_'$a1'2.fastq' -S $a1".sam"

gzip 'paired_trimmed_'$a1'1.fastq' 'paired_trimmed_'$a1'2.fastq'
singularity exec /usr/local/biotools/s/samtools:1.8--4 samtools view -bS $a1".sam" |singularity exec /usr/local/biotools/s/samtools:1.8--4 samtools sort -o $a1".bam"
singularity exec /usr/local/biotools/s/samtools:1.8--4 samtools index $a1".bam"
singularity exec /usr/local/biotools/b/bedtools:2.26.0gx--0 bedtools bamtobed -i $a1".bam" |sort -k 1,1 -k 2,2n > $a1".bed"
rm $a1".sam"
