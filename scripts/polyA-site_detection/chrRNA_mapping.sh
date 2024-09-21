#!/bin/sh
#$ -S /bin/sh
#$ -cwd

module load singularity
alias cutadapt='singularity exec /usr/local/biotools/c/cutadapt:1.15--py36_0 cutadapt'
alias STAR='singularity exec /usr/local/biotools/s/star:2.5.3a--0 STAR'
alias bedtools='singularity exec /usr/local/biotools/b/bedtools:2.26.0gx--0 bedtools'

for in_file in *_1.fastq.gz
do
file_name=$(echo $in_file| sed -e 's/.fastq.gz//')
#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o $file_name'_cutadapt.fastq.gz' $in_file  
#cutadapt -g "XT{151}" --minimum-length 40 --too-short-output $file_name'_too_short40_polyT151.fastq.gz' -o $file_name'_cutadapt_polyT151.fastq.gz' $file_name'_cutadapt.fastq.gz'
STAR --runThreadN 8 --genomeDir ~/refs/Arabidopsis_STAR_araport --readFilesCommand gunzip -c --readFilesIn $file_name'_cutadapt_polyT151.fastq.gz' --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFileNamePrefix ./$file_name --quantMode GeneCounts
bedtools genomecov -ibam $file_name'Aligned.sortedByCoord.out.bam' -g ~/refs/Chromosome_lengths_Arabidopsis.txt -bg -split -strand + > $file_name'_rev_strand.bedgraph'
bedtools genomecov -ibam $file_name'Aligned.sortedByCoord.out.bam' -g ~/refs/Chromosome_lengths_Arabidopsis.txt -bg -split -strand - > $file_name'_fwd_strand.bedgraph'
awk '{print $1,$2,$3,-$4}' $file_name'_rev_strand.bedgraph' > $file_name'_rev_strand1.bedgraph'
cat $file_name'_fwd_strand.bedgraph' $file_name'_rev_strand1.bedgraph' > $file_name'.bedgraph'
#awk '{print $1,$2,$3,$4/15.686927}' $file_name'.bedgraph' > $file_name'_unique_normalize.bedgraph' #uniquely mapped readsが15,686,927リードの場合
done

# unique map で normalize
awk '{print $1,$2,$3,$4/14.874329}' 20-702_1.bedgraph > 20-702_1_unique_normalize.bedgraph &
awk '{print $1,$2,$3,$4/21.749954}' 20-703_1.bedgraph > 20-703_1_unique_normalize.bedgraph &
awk '{print $1,$2,$3,$4/23.126465}' 20-704_1.bedgraph > 20-704_1_unique_normalize.bedgraph &
awk '{print $1,$2,$3,$4/24.054146}' 20-705_1.bedgraph > 20-705_1_unique_normalize.bedgraph &
awk '{print $1,$2,$3,$4/13.402709}' 20-706_1.bedgraph > 20-706_1_unique_normalize.bedgraph &
awk '{print $1,$2,$3,$4/17.955273}' 20-707_1.bedgraph > 20-707_1_unique_normalize.bedgraph &