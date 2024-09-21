## PEND! see also temp_chrRNA.R
# load plotlib function from https://github.com/shenlab-sinai/ngsplot/blob/develop/lib/plotlib.r
source('plotlib.r')

## ---function----
read_ngsplot_hmdata<-function(datapath){
  hm_list <- list()
  line_list <- list()
  agxtxt<-read.table(file.path(datapath,'avgprof.txt'),header=T,sep='\t')
  # Loop through each file
  for (i in 1:6) {
    # Read the file into the hm_list
    hm_list[[i]] <- read.table(file.path(datapath, paste0('hm', i, '.txt')), header=T, sep='\t')
    # Compute the mean and store it in line_list
    line_list[[i]] <- apply(hm_list[[i]][, 5:ncol(hm_list[[i]])], 2, mean, na.rm=T)
  }
  return(list("hm_list" = hm_list, "line_list" = line_list,"avgtxt" = agxtxt))
}

plotprofile<-function(line_list,columns_to_be_plotted=c(1,2,3,4),colors=c('gray70','lightblue2','gray40','dodgerblue3'),lwd=c(4,4,2,2),labelnames=colnames(agxtxt),abln=c(20,80)){
  ylims=c(min(unlist(line_list)),max(unlist(line_list)))
  for(i in 1:length(columns_to_be_plotted)){
    plot(1:101,line_list[[i]],type='l',col=colors[i],lwd=lwd[i],ylim=ylims,xlab='',ylab='',xaxt='n');par(new=T)
  }
  abline(v=abln)
  text(rep(75,4),seq(from=0.75*ylims[2],to=0.9*ylims[2],length.out=4),labelnames,col=colors,pos=4)
  
}


par(mar=c(0.2,2,0.2,0.2))

## remove when publish-----
system('mkdir ../data')
system('mkdir ../data/Figure6')
system('mkdir ../data/Figure6/atxr3_wt_chrRNA_body')
system('mkdir ../data/Figure6/atxr3_wt_rna_ChIP1_body')
system('mkdir ../data/Figure6/TSS_seq2_atxr3s__Rep2')
system('mkdir ../data/Figure6/TSS_seq2_atxr3s__Rep1')
system('mkdir ../data/Figure6/cdkc2_wt_H3K4me3')

system('cp /Users/Satoyo/data/chrRNA/ngsplots/atxr3_wt_chrRNA_body/*txt ../data/Figure6/atxr3_wt_chrRNA_body/')
system('cp /Users/Satoyo/data/ChIP1/RNA/atxr3_wt_rna_ChIP1_body/*txt ../data/Figure6/atxr3_wt_rna_ChIP1_body/')
system('cp /Users/Satoyo/data/TSS_seq2/metaplots/TSS_seq2_atxr3s_Rep2/*txt ../data/Figure6/TSS_seq2_atxr3s__Rep2/')
system('cp /Users/Satoyo/data/TSS_seq2/metaplots/TSS_seq2_atxr3s__Rep1/*txt ../data/Figure6/TSS_seq2_atxr3s__Rep1/')
system('cp /Users/Satoyo/data/chrRNA/ngsplots/cdkc2_wt_H3K4me3/*txt ../data/Figure6/cdkc2_wt_H3K4me3/')

##-----

# Fig.6a
# ------prep data --------
# avgprof.txt, hm*.txt files are the output of ngs.plot.r (https://github.com/shenlab-sinai/ngsplot/)
# ngs.plot.r was run with the following command
# the metaplots and heatmaps were generated with ngs.plot.r with the following commands
# ngs.plot.r -G Tair10 -R genebody -C atxr3_wt_allgene_bluegene_chrRNA.txt -O atxr3_wt_chrRNA_body -GO none -L 500 &
# ----atxr3_wt_allgene_bluegene_chrRNA.txt----
#20-702_1Aligned.sortedByCoord.out.bam /home/satoyo08/refs/gene_list/ChIP1_WT_expaverage_decending_order_nuclear.txt  "WT nuclear genes"
#20-707_1Aligned.sortedByCoord.out.bam /home/satoyo08/refs/gene_list/ChIP1_WT_expaverage_decending_order_nuclear.txt "atxr3 nuclear genes"
#20-702_1Aligned.sortedByCoord.out.bam /home/satoyo08/refs/gene_list/the_new_blue.txt  "WT ATXR3-marked"
#20-707_1Aligned.sortedByCoord.out.bam /home/satoyo08/refs/gene_list/the_new_blue.txt "atxr3 ATXR3-marked"
#20-702_1Aligned.sortedByCoord.out.bam /home/satoyo08/refs/gene_list/ChIP14_ATXR-tags_bound/ATXR3_TSS_AtID.txt "WT ATXR3-bound"
#20-707_1Aligned.sortedByCoord.out.bam /home/satoyo08/refs/gene_list/ChIP14_ATXR-tags_bound/ATXR3_TSS_AtID.txt "r3 ATXR3-bound"

#------visualize-------
datapath='../data/Figure6/atxr3_wt_chrRNA_body/'
DT<-read_ngsplot_hmdata(datapath)
plotprofile(DT$line_list,labelnames=c('WT all nuclear genes','atrx3 all nuclear genes','WT ATXR3-marked genes','atxr3 ATXR3-marked genes'))

# ----- Fig 6.b ------
#ngs.plot.r -G Tair10 -R genebody -C atxr3_wt_allgene_bluegene.txt -O atxr3_wt_rna_ChIP1_body -GO none -L 500 &
#merged_sorted_WT.bam /home/satoyo08/refs/gene_list/ChIP1_WT_expaverage_decending_order_nuclear.txt  "WT nuclear genes"
# merged_sorted_r3.bam  /home/satoyo08/refs/gene_list/ChIP1_WT_expaverage_decending_order_nuclear.txt "atxr3 nuclear genes"
# merged_sorted_WT.bam /home/satoyo08/refs/gene_list/the_new_blue.txt  "WT ATXR3-marked"
# merged_sorted_r3.bam  /home/satoyo08/refs/gene_list/the_new_blue.txt "atxr3 ATXR3-marked"
# merged_sorted_WT.bam /home/satoyo08/refs/gene_list/ChIP14_ATXR-tags_bound/ATXR3_TSS_AtID.txt "WT ATXR3-bound"
# merged_sorted_r3.bam /home/satoyo08/refs/gene_list/ChIP14_ATXR-tags_bound/ATXR3_TSS_AtID.txt "r3 ATXR3-bound"
datapath='../data/Figure6/atxr3_wt_rna_ChIP1_body'
DT<-read_ngsplot_hmdata(datapath)
plotprofile(DT$line_list,labelnames=c('WT all nuclear genes','atrx3 all nuclear genes','WT ATXR3-marked genes','atxr3 ATXR3-marked genes'))

#----- Fig. 6c -----
# ngs.plot.r -G Tair10 -R genebody -C cdkc2_wt_H3K4me3.txt -O cdkc2_wt_H3K4me3 -GO none -L 500 &
# /home/shusei89/ChIP_seq/Apr2022/fastp/371_K4me1_1.trim.bam /home/satoyo08/refs/gene_list/ChIP1_WT_expaverage_decending_order_nuclear.txt  "H3K4me1 WT nuclear genes"
# /home/shusei89/ChIP_seq/Apr2022/fastp/373_K4me1_1.trim.bam  /home/satoyo08/refs/gene_list/ChIP1_WT_expaverage_decending_order_nuclear.txt "cdkc2 nuclear genes"
# /home/shusei89/ChIP_seq/Apr2022/fastp/371_K4me1_1.trim.bam /home/satoyo08/refs/gene_list/the_new_blue.txt  "WT ATXR3-marked"
# /home/shusei89/ChIP_seq/Apr2022/fastp/373_K4me1_1.trim.bam  /home/satoyo08/refs/gene_list/the_new_blue.txt "cdkc2 ATXR3-marked"
# /home/shusei89/ChIP_seq/Apr2022/fastp/371_K4me1_1.trim.bam /home/satoyo08/refs/gene_list/ChIP14_ATXR-tags_bound/ATXR3_TSS_AtID.txt "WT ATXR3-bound"
# /home/shusei89/ChIP_seq/Apr2022/fastp/373_K4me1_1.trim.bam /home/satoyo08/refs/gene_list/ChIP14_ATXR-tags_bound/ATXR3_TSS_AtID.txt "cdkc2 ATXR3-bound"
datapath='../data/Figure6/cdkc2_wt_H3K4me3'
DT<-read_ngsplot_hmdata(datapath)
plotprofile(DT$line_list,labelnames=c('WT all nuclear genes','atrx3 all nuclear genes','WT ATXR3-marked genes','atxr3 ATXR3-marked genes'))



#---- Fig. 6d -----
# ngs.plot.r -G Tair10 -R tss -C config2.txt -O TSS_seq2_atxr3s__Rep1 -GO none -L 500 &
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_WT_4.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/ChIP1_WT_expaverage_decending_order_nuclear.txt  "WT_4 nuclear genes"
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_r3_2.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/ChIP1_WT_expaverage_decending_order_nuclear.txt "atxr3_2 nuclear genes"
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_WT_4.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/the_new_blue.txt  "WT_4 ATXR3-marked"
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_r3_2.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/the_new_blue.txt "atxr3_2 ATXR3-marked"
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_WT_4.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/ChIP14_ATXR-tags_bound/ATXR3_TSS_AtID.txt "WT_4 ATXR3-bound"
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_r3_2.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/ChIP14_ATXR-tags_bound/ATXR3_TSS_AtID.txt "atxr3_2 ATXR3-bound"
datapath='../data/Figure6/TSS_seq2_atxr3s__Rep1'
DT<-read_ngsplot_hmdata(datapath)
plotprofile(DT$line_list,labelnames=c('WT all nuclear genes','atrx3 all nuclear genes','WT ATXR3-marked genes','atxr3 ATXR3-marked genes'),abln=50)


# ---- Fig. 6e -----
# ngs.plot.r -G Tair10 -R tss -C config3.txt -O TSS_seq2_atxr3s_Rep2 -GO none -L 500 &
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_WT_2.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/ChIP1_WT_expaverage_decending_order_nuclear.txt  "WT_2 nuclear genes"
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_r3_3.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/ChIP1_WT_expaverage_decending_order_nuclear.txt "atxr3_3 nuclear genes"
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_WT_2.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/the_new_blue.txt  "WT_2 ATXR3-marked"
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_r3_3.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/the_new_blue.txt "atxr3_3 ATXR3-marked"
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_WT_2.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/ChIP14_ATXR-tags_bound/ATXR3_TSS_AtID.txt "WT_2 ATXR3-bound"
# Dr_Inagaki_TSS_rmdup_BAM/Dr_Inagaki_r3_3.UMI_1mm_deduplicated_sorted.bam /home/satoyo08/refs/gene_list/ChIP14_ATXR-tags_bound/ATXR3_TSS_AtID.txt "atxr3_3 ATXR3-bound"
datapath='../data/Figure6/TSS_seq2_atxr3s__Rep2'
DT<-read_ngsplot_hmdata(datapath)
plotprofile(DT$line_list,labelnames=c('WT all nuclear genes','atrx3 all nuclear genes','WT ATXR3-marked genes','atxr3 ATXR3-marked genes'),abln=50)

