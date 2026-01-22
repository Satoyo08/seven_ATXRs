# seven_ATXRs

Datasets and codes for "Comparative characterization of chromatin-targeting mechanisms across seven H3K4 methyltransferases in Arabidopsis" in Plant Communication (2026)

# codes for data preprocessing (in the 'scripts' folder)

## ChIP-seq

For ATX(R)-FLAG datasets;

For both of the two replicates, raw data files (fastq.gz) were trimmed by `scripts/trimmer_array.sh`, then mapped and counted by `scripts/fastq_to_bed_array.sh`. The resulting coverage files (bed) were then analyzed in R scripts in the folder the_seven_ATXRs.

For H3K4me ChIP-seq of mutants of *atx3, atx4, atx5* (Fig. 4);

Used read1 of the pair end reads. For both of the two replicates, raw data files (fastq.gz) were processed using `scripts/modified_from_SI_bowtie_loop_all_random_ChIP_no_singularity.sh`

Raw data files and browser track files (bw) are provided at GSE277439

## Splicing Analysis

You can find the outline in the supplemental material of the paper.
For the scripts and the key intermediate data files, see /scripts/vast-tools and /data/vast-tools. 

## polyA-site detection

You can find the outline in the supplemental material of the paper.
For the scripts and the key intermediate data files, see /scripts/polyA-site_detection for the scripts and /data/polyA-site_detection.


## AlphaScreen 

1490 Arabidopsis TFs were tested against ATX3 protein […] The full result is provided in Supplementary Table N in the paper.

The binding sequence of screened 1490 TFs was extracted from planttfdb (Ath_TF_binding_motifs.meme).  The genomic motif enrichment of all TFBS in the database was calculated with ame, comparing the ATX3-bound TSS regions to ATX3-unbound TSS regions. The distribution of the AlphaScreen Score and enrichment score (p-val) are provided in Supplementary Figure 8. The TFs that passed either of the following two criteria were presented in Figure 7: AlphaScreen score of ATX3/WGE (background control) > 7, or ame motif enrichment score -log10(pval) > 230. 


## TF search of SVM motif

Notable motifs discovered through SVM modeling were screened for the matching TFBS in the database. First, the query sequences, ‘RGCCCAW,’ 'TCGTCGTC,' 'TCYGATTC', and 'SCGGCGR' were converted to position matrix in meme format by `scripts/make_TOMTOM_query.ipyenv`. Derived meme files were searched against the following databases using tomotom: JASPAR_CORE_2014_plants.meme, ArabidopsisDAPv1.meme and ArabidopsisPBM_20140210.meme in the MEME Suite’s motif_databases.12.19, and Ath_TF_binding_motifs.meme in PlantTFDB.  

# Figure-by-figure description of data visualization

### Figure 1

AlphaFold models were obtained from UniProt on May 18, 2024, and are available at data/alpha_fold/pbd. Domain annotations file can be found in the GitHub repo for our previous paper (https://github.com/Satoyo08/Arabidopsis_H3K4me1/tree/main/data/Figure1/*csv), which were also obtained from UniProt.
Models were colored and visualized using /scripts/chimera_commands.txt


The heatmaps were generated using deeptools. For genomic localizations of H3K4-ATX(R)s...
```
deeptools bamCompare -b1 $sample -b2 $control -o $outname".bw"  -p 8
deeptools computeMatrix scale-regions -p 8 -S $bwfiles -R araport11_all_sorted_by_RNA_RPKM.bed -b 500 -a 500 -o $bwname"_over_"$roiname".mat.gz"
deeptools plotHeatmap -p --sortRegions no --heatmapHeight 14 -m $bwname"_over_"$roiname".mat.gz" -o $bwname"_over_"$roiname".png" 
```
For H3K4me profiles, the ChIP-seq data was reanalyzed from our previous paper (Oya et al., 2022)

```
#WT and mutant profile
deeptools bamCoverage -b $bamfile -o $bwname".bw"
deeptools computeMatrix scale-regions -p 8 -S $bwfiles -R araport11_all_sorted_by_RNA_RPKM.bed -b 500 -a 500 -o $bwname"_over_"$roiname".mat.gz"
deeptools plotHeatmap -p 8 --sortRegions no --heatmapHeight 14 -m $bwname"_over_"$roiname".mat.gz" -o $bwname"_over_"$roiname".png" 

# mutant/WT profile
deeptools bamCompare -b1 $sample -b2 $control -o $outname".bw"  -p 8
deeptools computeMatrix scale-regions -p 8 -S $bwfiles -R araport11_all_sorted_by_RNA_RPKM.bed -b 500 -a 500 -o $bwname"_over_"$roiname".mat.gz"
deeptools plotHeatmap  --zMin -0.8 --zMax 0.8 --sortRegions no  --colorMap YlGnBu --heatmapHeight 14 -m $mat.gz -o $png
```

### Figure 2

Full codes for training and visualization with links to the models and data can be found in Rscripts/Figure2.r

### Figure 3

`scripts/SVM_training.ipyenv` was used to train SVM models.
Full codes for visualization and the links to the models can be found in Rscripts/Figure3.R

### Figure 4

ChIP-seq peaks were called by ` macs2 callpeak -t mutant.bam -c control.bam -n output file -f BAM -g 1.3e8`. For a-d, full codes for visualization and the links to the intermediate files can be found in Figure4_SFig_5_6.R

### Figure 5 

See Figure/table legend, material and method section and Supplementary tables 2 and 3 in the paper offers enough explanation.


### Figure 6

See the Figure/table legend, material and method section and Supplementary tables 4 in the paper.

### Supplementary Figure 1

The heatmaps were generated using deeptools as described in the above section for Figure 1. Codes for venn diagram and hypergeometric test are provided in Supplementary_Figure1.R

### Supplementary Figure 2
Metaplots were generated using ngs.plot.r, and re-formatted in R by log-transforming the y-axis. The intermediate files (output datafile from ngs.plot.r) for metaplots can be found at data/RF_feature_metaplot

### Supplementary Figure 3
Full codes for training and visualization with links to the models and data can be found in Rscripts/Figure2.r

### Supplementary Figure 4,5,6 

Full codes for visualization and the links to the models can be found in Rscripts/Figure3.R

### Supplementary Figure 7,8

Full codes for visualization and the links to the intermediate files can be found in Rscripts/Figure4_SFig_7_8.R

### Supplementary Figure 10

Full codes for visualization and the links to the intermediate files can be found in Rscripts/Supplementary_Figure10.r

### Supplementary Figure 12

See the Supplementary Methods and Splicing Analysis section above.

### Supplementary Figure 13

See the Supplementary Methods and polyA-site detection section above.
For b, Full codes for visualization and the links to the intermediate files can be found in Rscripts/Supplementary_Figure13.r

### Supplementary Figure 15

Full codes for visualization and the links to the intermediate files can be found in Rscripts/ATX345_GO.r
