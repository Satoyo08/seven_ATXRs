## Align

First, the pair-end RNA-seq data was aligned with the following script;

```
#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd

# vast_tools_align.sh
fastq_array=(`cut -f 1 $1`)
pair1=${fastq_array[$SGE_TASK_ID]}
pair2=$(echo $pair1 | sed -e 's/_1\./_2\./')

alias vast-tools='~/bin/vast-tools/vast-tools'
vast-tools align $pair1 $pair2 -sp Ath --bowtieProg /home/satoyo08/bin/bowtie-1.0.0/bowtie
```

`qsub -t 1-84:1 vast_tools_align.sh fastq_list_R1 (-j 17116848) # align`

## Merge

RNA-seq data consists of three biological replicates (RNA from different plant individuals) for WT and atxr7 mutants. Aligned data was merged for each genotype for further analysis, as done so in previous research (Martin *et al.*, 2021). 

```
module load r
vast-tools merge --sp araTha10 --groups v1_group_config 
```

```
### v1_group_config
# 109 = atxr7, 122 = WT. 109_1,_2,_3 stands for three biological replicates.
v1_109_1_1	v1_109
v1_109_2_1	v1_109
v1_109_3_1	v1_109
v1_122_1_1	v1_122
v1_122_2_1	v1_122
v1_122_3_1	v1_122
```



## Detection of differentially spliced sites

```
V1_merged_tab=INCLUSION_LEVELS_FULL-araTha10-7.tab
vast-tools tidy $V1_merged_tab --min_N 1

vast-tools compare $V1_merged_tab -a v1_122 -b v1_109 --min_dPSI 15  --min_range 5 --GO --sp Ath
vast-tools compare $V1_merged_tab -a v1_109 -b v1_122 --min_dPSI 15  --min_range 5 --GO --sp Ath
```

These process generated two files,  `DiffAS-araTha10-7-dPSI15-range5-min_ALT_use25-upreg_ALT_v1_122-vs-v1_109.tab` and `DiffAS-araTha10-7-dPSI15-range5-min_ALT_use25-upreg_ALT_v1_109-vs-v1_122.tab`. The former file represents an alternative splicing event that increased in the atxr7 mutant and the latter opposite.



## supplementarily data

The following intermediate files are available in the same folder as this script.

`DiffAS-araTha10-7-dPSI15-range5-min_ALT_use25-upreg_ALT_v1_122-vs-v1_109.tab` 

 `DiffAS-araTha10-7-dPSI15-range5-min_ALT_use25-upreg_ALT_v1_109-vs-v1_122.tab`

`INCLUSION_LEVELS_FULL-araTha10-7.tab`

