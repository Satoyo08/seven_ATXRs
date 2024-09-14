# ChIP-seq

For ATX(R)-FLAG datasets;

For both of the two replicates, raw data files (fastq.gz) were trimmed by `trimmer_array.sh`, then mapped and counted by  `fastq_to_bed_array.sh`. The resulting coverage files (bed) were then analyzed in R scripts in the folder the_seven_ATXRs.

For H3K4me ChIP-seq of mutants of *atx3, atx4, atx5* (Fig. 4)

Used read1 of the pair end reads. For both of the two replicates, raw data files (fastq.gz) were processed using `modified_from_SI_bowtie_loop_all_random_ChIP_no_singularity.sh`

Raw data files and browser track files (bw) are provided at [GEOaccession]

## splicing analysis

see /scripts/vast-tools for the key intermediate data and the scripts

## polyA-site detection

see /scripts/polyA-site_detection for the scripts and /data/polyA for the key intermediate data.





# AlphaScreen and it’s analysis

1490 Arabidopsis TFs were tested against ATX3 protein […] The full result is provided in Supplementary Table N. 

The binding sequence of screened 1490 TFs were extracted from planttfdb (Ath_TF_binding_motifs.meme).  The genomic motif enrichment of all TFBS in the database was calculated with ame, comparing the ATX3-bound TSS regions to ATX3-unbound TSS regions. The distribution of the AlphaScreen Score and enrichment score (p-val) are provided in Supplementary Figure 8. The TFs that passed either of the following two criteria were presented in Figure 7; AlphaScreen score of ATX3/WGE (background control) > 7 , or ame motif enrichment score -log10(pval) > 230. 



# TF search of SVM motif

Notable motifs discovered through SVM modeling were screened for the matching TFBS in the database. First, the query sequences, ‘RGCCCAW’,'TCGTCGTC','TCYGATTC' and 'SCGGCGR' were converted to position matrix in meme format by `make_TOMTOM_query.ipyenv` . Derived meme files were searched against the following databases using tomotom; JASPAR_CORE_2014_plants.meme, ArabidopsisDAPv1.meme and ArabidopsisPBM_20140210.meme in the MEME Suite’s motif_databases.12.19, and Ath_TF_binding_motifs.meme in PlantTFDB. 



# Fig_by_Fig scripts

##### Fig1

- [ ] AlphaFold; not necessary? but available anyway.
- [ ] Deeptools
- [ ] Deeptools-difference

##### Fig2

Done. script and models and the data frame.

##### Fig3

- [x] ROC. String for rep1,2 
- [x] Heatmap
- [x] SVM scripts
  - [x] ipyenv in  https://github.com/Satoyo08/Arabidopsis_H3K4me1 was used. Edited and moved to the shells/


##### Fig4

- [x] barplot peaks, scatter plot, box plot -> Figure4_SFig_5_6.R
- [ ] Western blot
- [ ] 



