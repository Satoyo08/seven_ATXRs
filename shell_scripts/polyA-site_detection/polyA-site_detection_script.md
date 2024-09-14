# PolyA analysis

### Dataset1

The sequences were adaptor-trimmed, polyA-extracted and mapped by `array_cutadapt_polyT_STAR.sh`

The reads were cleaned up by removing duplicates and compared with the annotated polyA sites downloaded from APAdb(`arabidopsis_thaliana.high_confidence.PAC_bed6.bed`)

```
polyA_bed=arabidopsis_thaliana.high_confidence.PAC_bed6.bed

for inbam in trimmed_polyT_*Aligned.sortedByCoord.out.bam
do
prefix=$(echo $inbam | sed -e 's/Aligned\.sortedByCoord\.out\.bam/_/')
samtools sort $inbam -o positionsort.bam
samtools markdup -r -s positionsort.bam $prefix"_rmdup.bam"
bedtools bamtobed -i $prefix"_rmdup.bam" > $prefix".bed"
bedtools flank -s -i $prefix".bed" -g $glen -l 1 -r 0 |sort -k 1,1 -k 2,2n> $prefix"5end.bed"
bedtools closest -S -D a -a $prefix"5end.bed" -b $polyA_bed > $prefix"5end_intersect_polyA.bed"
done
```

The $prefix"_rmdup.bam" files were pooled based on genotype by samtools merge to generate polyT_WT.bam and polyT_ATXR7.bam. Metaplot was generated from the merged bamfiles were plotted by the following ngs.plot.r command;

`ngs.plot.r -G Tair10 -R genebody -C ngs_polyA_config.txt -O ATXR7_WT_polyA -GO none -L 500`

```
polyT_WT.bam ChIP1_WT_expaverage_decending_order_nuclear.txt "WT"
polyT_ATXR7.bam ChIP1_WT_expaverage_decending_order_nuclear.txt "atxr7"
polyT_WT.bam top_3000_genelist_ATXR7_ChIP8_TESAtID_sorted_bound.bed.txt "WT_ATXR7-bound"
polyT_ATXR7.bam top_3000_genelist_ATXR7_ChIP8_TESAtID_sorted_bound.bed.txt "atxr7_ATXR7_bound"
```

### Dataset2

The adaptors were trimmed by `array_cutadapt_rm_adapter.sh`. Reads with polyA were extracted by `array_cutadapt_polyT.sh`. Then the extracted reads were mapped by `star_polyT.sh`.

The reads were cleaned up by removing duplicates and compared with the annotated polyA sites downloaded from APAdb (`arabidopsis_thaliana.high_confidence.PAC_bed6.bed`)

```
polyA_bed=arabidopsis_thaliana.high_confidence.PAC_bed6.bed

for inbam in trimmed_polyT_*Aligned.sortedByCoord.out.bam
do
prefix=$(echo $inbam | sed -e 's/Aligned\.sortedByCoord\.out\.bam/_/')
samtools collate -o namecollate.bam $inbam 
samtools fixmate -m namecollate.bam fixmate.bam
samtools sort -o positionsort.bam fixmate.bam
samtools markdup -r -s positionsort.bam markdup.bam
samtools view -bf 0x40 markdup.bam > $prefix"R1.bam" 
bedtools bamtobed -i $prefix"R1.bam" > $prefix"R1.bed"
bedtools flank -s -i $prefix"R1.bed" -g $glen -l 1 -r 0 |sort -k 1,1 -k 2,2n> $prefix"R1_5end.bed"
bedtools closest -S -D a -a $prefix"R1_5end.bed" -b $polyA_bed > $prefix"R1_5end_intersect_polyA.bed"
done
```

the $prefix"R1.bam" files were pooled based on genotype by samtools merge to generate R1_of_WT.bam and R1_of_ATXR7.bam. Metaplot was generated from the merged bamfiles were plotted by the following ngs.plot.r command;

`ngs.plot.r -G Tair10 -R genebody -C ngs_polyA_config.txt -O ATXR7_WT_polyA -GO none -L 500`

```
#ngs_polyA_config.txt
R1_of_WT.bam ChIP1_WT_expaverage_decending_order_nuclear.txt "WT"
R1_of_ATXR7.bam ChIP1_WT_expaverage_decending_order_nuclear.txt "atxr7"
R1_of_WT.bam top_3000_genelist_ATXR7_ChIP8_TESAtID_sorted_bound.bed.txt "WT_ATXR7-bound"
R1_of_ATXR7.bam top_3000_genelist_ATXR7_ChIP8_TESAtID_sorted_bound.bed.txt "atxr7_ATXR7_bound"
```





### reference table of APA db and gene annotation

```
genes=araport11_all_sorted_nuclear.bed
bedtools intersect -wao -loj -s -a $polyA_bed -b $genes > araport_genes_and_polyA.bed
grep "AT" ../araport_genes_and_polyA.bed | sort -k 11,11n -k 2,2n > ../araport_genes_and_polyA_gonly.bed
```

Downstream analysis was conducted in R (polyA_stats.R)

