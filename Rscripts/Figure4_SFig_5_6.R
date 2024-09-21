source('custom_functions.r')

### barplot Number of peaks ---
# (Rep2) read peak files (MACS2 output narrowPeak files)
atx3_me3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep2/K4me3_atx3-1_reduction_peaks.narrowPeak",header=F)
atx45_me3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep2/K4me3_atx45_reduction_peaks.narrowPeak",header=F)
atx345_me3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep2/K4me3_atx3-145_reduction_peaks.narrowPeak",header=F)
atx3_me2_R<-read.table("../data/Figure4_SFig5_6/peaks_rep2/K4me2_atx3-1_reduction_peaks.narrowPeak",header=F)
atx45_me2_R<-read.table("../data/Figure4_SFig5_6/peaks_rep2/K4me2_atx45_reduction_peaks.narrowPeak",header=F)
atx345_me2_R<-read.table("../data/Figure4_SFig5_6/peaks_rep2/K4me2_atx3-145_reduction_peaks.narrowPeak",header=F)
atx3_H3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep2/H3_atx3-1_reduction_peaks.narrowPeak",header=F)#empty
atx3_H3_R<-data.frame(0);colnames(atx3_H3_R)<-'V7'
atx45_H3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep2/H3_atx45_reduction_peaks.narrowPeak",header=F)
atx345_H3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep2/H3_atx3-145_reduction_peaks.narrowPeak",header=F)
LIST1<-list(atx3_me3_R,atx45_me3_R,atx345_me3_R,atx3_me2_R,atx45_me2_R,atx345_me2_R,atx3_H3_R,atx45_H3_R,atx345_H3_R)

# (Rep1) read peak files (MACS2 output narrowPeak files)
atx3_me3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep1/H3K4me3_274_reduction_peaks.narrowPeak",header=F)
atx45_me3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep1/H3K4me3_273_reduction_peaks.narrowPeak",header=F)
atx345_me3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep1/H3K4me3_275_reduction_peaks.narrowPeak",header=F)
atx3_me2_R<-read.table("../data/Figure4_SFig5_6/peaks_rep1/H3K4me2_274_reduction_peaks.narrowPeak",header=F)
atx45_me2_R<-read.table("../data/Figure4_SFig5_6/peaks_rep1/H3K4me2_273_reduction_peaks.narrowPeak",header=F)
atx345_me2_R<-read.table("../data/Figure4_SFig5_6/peaks_rep1/H3K4me2_275_reduction_peaks.narrowPeak",header=F)
atx3_H3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep1/H3_274_reduction_peaks.narrowPeak",header=F)
atx45_H3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep1/H3_273_reduction_peaks.narrowPeak",header=F)
atx345_H3_R<-read.table("../data/Figure4_SFig5_6/peaks_rep1/H3_275_reduction_peaks.narrowPeak",header=F) # empty
atx345_H3_R<-data.frame(0);colnames(atx345_H3_R)<-'V7'
LIST2<-list(atx3_me3_R,atx45_me3_R,atx345_me3_R,atx3_me2_R,atx45_me2_R,atx345_me2_R,atx3_H3_R,atx45_H3_R,atx345_H3_R)


LIST<-LIST1
selected_peaks<-c()
cutoff=2.5
for(i in 1:9){
  DF<-LIST[[i]]
  selected_peaks[i]<-length(DF[DF$V7>cutoff,1])
}
peaks<-data.frame(name=c("x3_me3","x45_me3","x345_me3","x3_me2","x45_me2","x345_me2","x3_H3","x45_H3","x345_H3"));peaks$rep2<-selected_peaks

LIST<-LIST2
selected_peaks<-c()
for(i in 1:9){
  DF<-LIST[[i]]
  selected_peaks[i]<-length(DF[DF$V7>cutoff,1])
}
peaks$rep1<-selected_peaks

## plot barplot 
peaks$mean<-apply(peaks[,2:3],1,mean)
PM<-as.matrix(cbind(peaks$mean[7:9],peaks$mean[4:6],peaks$mean[1:3]))
barplot(t(PM),beside=TRUE,ylim=c(0,8400))
x<-c(1:3,5:7,9:11)+0.5
points(x,peaks$rep1[c(7,4,1,8,5,2,9,6,3)],pch=16)
points(x,peaks$rep2[c(7,4,1,8,5,2,9,6,3)],pch=1)



### ---- (Rep1) scatter plot 3x3 ------
# read ChIP-seq files
ChIP16_RPM<-read.table("../data/Figure4_SFig5_6/ChIP16_RPKM.txt",header=TRUE)
ChIP16_RPKM<-read.table("../data/Figure4_SFig5_6/ChIP16_RPKM.txt",header=TRUE)

# read gene list--
dir2='../data/refs/ATX(R)s_bound_gens_list/rep1'
files<-list.files(dir2,full.names = TRUE);files
ATX3_bound<-read.table(files[3])$V1
ATX4_bound<-read.table(files[4])$V1
ATX5_bound<-read.table(files[5])$V1

AL<-merge(ChIP16_RPKM,ChIP16_RPM,by='ID');AL[,ncol(AL)+1]<-0
# -- not used for plotting -- 
AL[AL$ID %in% new_red,]$V26<-AL[AL$ID %in% new_red,]$V26+1
AL[AL$ID %in% new_green,]$V26<-AL[AL$ID %in% new_green,]$V26+2
AL[AL$ID %in% new_blue,]$V26<-AL[AL$ID %in% new_blue,]$V26+4
# ----

AL[,ncol(AL)+1]<-0
AL[AL$ID %in% ATX3_bound,]$V27<-AL[AL$ID %in% ATX3_bound,]$V27+1
AL[AL$ID %in% ATX4_bound,]$V27<-AL[AL$ID %in% ATX4_bound,]$V27+2
AL[AL$ID %in% ATX5_bound,]$V27<-AL[AL$ID %in% ATX5_bound,]$V27+4

pgy_pal<-c("gray80",pallets[c(3,10,6,11,7)],"#2D004B","black")

wt_col<-c(2,2,2,18,18,18,22,22,22);colnames(AL)[wt_col] # H3 for RPKM, H3K4me2, me3 for RPM. 
mutant_col<-c(4,3,5,20,19,21,24,23,25);colnames(AL)[mutant_col]

dev.off();par(mfrow=c(3,3),mar=c(4,4,2,2),cex.lab=1.2)
for(i in 1:length(wt_col)){
  plot_shu(AL[,wt_col[i]],AL[,mutant_col[i]],col=pgy_pal[AL[,ncol(AL)] + 1],
           #main=colnames(AL)[mutant_col[i]],
           pch=20,lwd=0.01) 
}
# venn diagram
table(AL$V27)
library(eulerr)
fit1 <- euler(c(
  "A" = 1679,
  "B" = 1312,
  "C" = 1532,
  "A&B" = 568,
  "A&C" = 348,
  "B&C" = 715,
  "A&B&C" = 405
))
plot(fit1) #adjust color in affinity designer


# 1. boxplot
dev.off()
AL$K4me3_atx345_div_atx45<-AL$K4me3_atx3.145.y-AL$K4me3_atx45.y # RPM/RPM
cpl<-c('gray80',pgy_pal[c(2,3,5)])
ALg<-gonly(AL)
boxplot(ALg$K4me3_atx345_div_atx45,
        AL[AL$ID %in% ATX3_bound,]$K4me3_atx345_div_atx45,
        AL[AL$ID %in% ATX4_bound,]$K4me3_atx345_div_atx45,
        AL[AL$ID %in% ATX5_bound,]$K4me3_atx345_div_atx45,
        notch=TRUE,outline=F,names=c('protein_coding','ATX3_bound','ATX4_bound','ATX5_bound'),las=2,col=cpl
        ,ylab='AL$K4me3_atx3.145.y-AL$K4me3_atx45.y')



### (Rep2) scatter plot 3x3 ---------
ChIP17_RPM<-read.table("../data/Figure4_SFig5_6/ChIP17_RPM.txt",header=TRUE)
ChIP17_RPKM<-read.table("../data/Figure4_SFig5_6/ChIP17_RPKM.txt",header=TRUE)


AL<-merge(ChIP17_RPKM,ChIP17_RPM,by='ID');AL[,ncol(AL)+1]<-0
# -- not used for plotting -- 
AL[AL$ID %in% new_red,]$V26<-AL[AL$ID %in% new_red,]$V26+1
AL[AL$ID %in% new_green,]$V26<-AL[AL$ID %in% new_green,]$V26+2
AL[AL$ID %in% new_blue,]$V26<-AL[AL$ID %in% new_blue,]$V26+4
# ---- 

AL[,ncol(AL)+1]<-0
AL[AL$ID %in% ATX3_bound,]$V27<-AL[AL$ID %in% ATX3_bound,]$V27+1
AL[AL$ID %in% ATX4_bound,]$V27<-AL[AL$ID %in% ATX4_bound,]$V27+2
AL[AL$ID %in% ATX5_bound,]$V27<-AL[AL$ID %in% ATX5_bound,]$V27+4


wt_col<-c(2,2,2,18,18,18,22,22,22);colnames(AL)[wt_col] # H3 for RPKM, H3K4me2, me3 for RPM. 
mutant_col<-c(4,3,5,20,19,21,24,23,25);colnames(AL)[mutant_col]

dev.off();par(mfrow=c(3,3),mar=c(4,4,2,2),cex.lab=1.2)
for(i in 1:length(wt_col)){
  plot_shu(AL[,wt_col[i]],AL[,mutant_col[i]],col=pgy_pal[AL[,ncol(AL)] + 1],
           #main=colnames(AL)[mutant_col[i]],
           pch=20,lwd=0.01) 
}


#  (rep2) boxplot
dev.off()
AL$K4me3_atx345_div_atx45<-AL$K4me3_atx3.145.y-AL$K4me3_atx45.y # RPM/RPM
cpl<-c('gray80',pgy_pal[c(2,3,5)])
ALg<-gonly(AL)
boxplot(ALg$K4me3_atx345_div_atx45,
        AL[AL$ID %in% ATX3_bound,]$K4me3_atx345_div_atx45,
        AL[AL$ID %in% ATX4_bound,]$K4me3_atx345_div_atx45,
        AL[AL$ID %in% ATX5_bound,]$K4me3_atx345_div_atx45,
        notch=TRUE,outline=F,names=c('protein_coding','ATX3_bound','ATX4_bound','ATX5_bound'),las=2,col=cpl
        ,ylab='AL$K4me3_atx3.145.y-AL$K4me3_atx45.y')



