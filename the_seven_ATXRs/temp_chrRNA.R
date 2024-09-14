setwd("/Users/Satoyo/data/chrRNA/")
hmplist<-list.files("/Users/Satoyo/data/chrRNA/rep2/atxr3_wt_chrRNA_body_Rep2",pattern = "hm",full.names = TRUE)

dev.off()
colpal<-c("gray70","skyblue","none","none","gray40","navy")
for(i in c(1,2,5,6)){
  hmp1<-read.table(hmplist[i],fill=T,header=T)
  summy<-apply(log(hmp1[,5:ncol(hmp1)]+0.01),2,mean)
  plot(1:length(summy),summy,type='l',col=colpal[i],lwd=2,ylim=c(-2,1));par(new=T)
}

legend("bottom",legend=c("WT nuclear coding","atxr3 nuclear coding","WT ATXR3-marked","atxr3 ATXR3-marked"),fill=colpal[c(1,2,5,6)])

### bedgraph 集計

#rep1
WT_rep1<-read.table("/Users/Satoyo/data/chrRNA/bedgraph/sense_20-702_1_unique_normalize.bedgraph")
atxr3_rep1<-read.table("/Users/Satoyo/data/chrRNA/bedgraph/sense_20-707_1_unique_normalize.bedgraph")
WT_n<-data.frame(abs(as.numeric(as.character(WT_rep1$V8))),"WT_nuc");colnames(WT_n)<-c("Value","Category")
r3_n<-data.frame(abs(as.numeric(as.character(atxr3_rep1$V8))),"r3_nuc");colnames(r3_n)<-c("Value","Category")
WT_marked<-data.frame(abs(as.numeric(as.character(WT_rep1[WT_rep1$V4 %in% new_blue,]$V8))),"WT_ATXR3_marked");colnames(WT_marked)<-c("Value","Category")
r3_marked<-data.frame(abs(as.numeric(as.character(atxr3_rep1[atxr3_rep1$V4 %in% new_blue,]$V8))),"r3_ATXR3_marked");colnames(r3_marked)<-c("Value","Category")
data<-rbind(rbind(rbind(WT_n,r3_n),WT_marked),r3_marked)

library(ggplot2)

# Generate the violin plot
ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_violin(trim = FALSE) +
  scale_y_log10()+
  labs(title = "Rep2 log10 sense", x = "Category", y = "Value") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black")+
  theme_minimal()

# rep2
WT_rep1<-read.table("/Users/Satoyo/data/chrRNA/rep2/sense_chrRNA_352_1_unique_normalize.bedgraph")
atxr3_rep1<-read.table("/Users/Satoyo/data/chrRNA/rep2/sense_chrRNA_353_1_unique_normalize.bedgraph")
WT_n<-data.frame(abs(as.numeric(as.character(WT_rep1$V7))),"WT_nuc");colnames(WT_n)<-c("Value","Category")
r3_n<-data.frame(abs(as.numeric(as.character(atxr3_rep1$V7))),"r3_nuc");colnames(r3_n)<-c("Value","Category")
WT_marked<-data.frame(abs(as.numeric(as.character(WT_rep1[WT_rep1$V4 %in% new_blue,]$V7))),"WT_ATXR3_marked");colnames(WT_marked)<-c("Value","Category")
r3_marked<-data.frame(abs(as.numeric(as.character(atxr3_rep1[atxr3_rep1$V4 %in% new_blue,]$V7))),"r3_ATXR3_marked");colnames(r3_marked)<-c("Value","Category")
data<-rbind(rbind(rbind(WT_n,r3_n),WT_marked),r3_marked)

## mRNAをlogに(see also pTEFb_and_atxr3.md)
dir<-'/Users/Satoyo/data/ChIP1/RNA/atxr3_wt_rna_ChIP1_body'
hmplist<-list.files(dir,pattern = "hm",full.names = TRUE)
dev.off()
colpal<-c("gray70","skyblue","gray40","navy","none","none")
for(i in c(1,2,3,4)){
  hmp1<-read.table(hmplist[i],fill=T,header=T)
  summy<-apply(log(hmp1[,5:ncol(hmp1)]+0.01),2,mean)
  plot(1:length(summy),summy,type='l',col=colpal[i],lwd=2,ylim=c(-4.7,1.5),ylab='log(+0.01)');par(new=T)
}
legend("bottom",legend=c("WT nuclear coding","atxr3 nuclear coding","WT ATXR3-marked","atxr3 ATXR3-marked"),fill=colpal[c(1,2,3,4)])

## TSS-seq をlogに(TSS_Seq_for_ATXR3.md)
dir<-'/Users/Satoyo/data/TSS_seq2/metaplots/TSS_seq2_atxr3s__Rep1'
hmplist<-list.files(dir,pattern = "hm",full.names = TRUE)
dev.off()
colpal<-c("gray70","skyblue","gray40","navy","none","none")
for(i in c(1,2,3,4)){
  hmp1<-read.table(hmplist[i],fill=T,header=T)
  summy<-apply(log(hmp1[,5:ncol(hmp1)]+0.01),2,mean)
  plot(1:length(summy),summy,type='l',col=colpal[i],lwd=2,ylim=c(-4.7,4),ylab='log(+0.01)');par(new=T)
}
legend("topleft",legend=c("WT nuclear coding","atxr3 nuclear coding","WT ATXR3-marked","atxr3 ATXR3-marked"),fill=colpal[c(1,2,3,4)])

dev.off()
colpal<-c("gray70","skyblue","gray40","navy","none","none")
for(i in c(1,2,3,4)){
  hmp1<-read.table(hmplist[i],fill=T,header=T)
  summy<-apply(hmp1[,5:ncol(hmp1)],2,mean)
  plot(1:length(summy),summy,type='l',col=colpal[i],lwd=2,ylim=c(0,120),ylab='',xlab='');par(new=T)
}
legend("topleft",legend=c("WT nuclear coding","atxr3 nuclear coding","WT ATXR3-marked","atxr3 ATXR3-marked"),fill=colpal[c(1,2,3,4)])





dir<-'/Users/Satoyo/data/TSS_seq2/metaplots/TSS_seq2_atxr3s_Rep2'
hmplist<-list.files(dir,pattern = "hm",full.names = TRUE)
dev.off()
colpal<-c("gray70","skyblue","gray40","navy","none","none")
for(i in c(1,2,3,4)){
  hmp1<-read.table(hmplist[i],fill=T,header=T)
  summy<-apply(log(hmp1[,5:ncol(hmp1)]+0.01),2,mean)
  plot(1:length(summy),summy,type='l',col=colpal[i],lwd=2,ylim=c(-4.7,4),ylab='log(+0.01)');par(new=T)
}
legend("bottom",legend=c("WT nuclear coding","atxr3 nuclear coding","WT ATXR3-marked","atxr3 ATXR3-marked"),fill=colpal[c(1,2,3,4)])

# TSS, mRNA and chrRNA results are funky in marked or bound genes.
# how about... re-doing in marked and bound genes?

path_rep1<-"/Users/Satoyo/Desktop/sevenATXRpaper/seven_ATXRs_localization_git/data/refs/ATX(R)s_bound_gens_list/rep1/ATXR3_TSS_AtID.txt"
ATXR3_bound<-read.table(path_rep1,header=F)
bound_marked<-intersect(new_blue,ATXR3_bound$V1)
length(bound_marked)
write.table(bound_marked,file="/Users/Satoyo/data/chrRNA/bound_marked.txt",quote=F,row.names = F,col.names = F)
# -> ngs plot, see chrRNA.md. 全部 /data/chrRNA に DL

setwd("/Users/Satoyo/data/chrRNA/ngsplots/")
dir<-'atxr3_wt_rna_ChIP1_body/'
hmplist<-list.files(dir,pattern = "hm",full.names = TRUE)
dev.off()
par(mfrow=c(2,1))
colpal<-c("gray70","skyblue","gray40","navy")
ymax=37
for(i in c(1,2,3,4)){
  hmp1<-read.table(hmplist[i],fill=T,header=T)
  summy<-apply(hmp1[,5:ncol(hmp1)],2,mean)
  plot(1:length(summy),summy,type='l',col=colpal[i],lwd=2,ylim=c(0,ymax),ylab='read count',xlab="");par(new=T)
}
mtext("mRNA")
legend("topleft",legend=c("WT nuc","atxr3 nuc","WT ATXR3-mb","atxr3 ATXR3-mb"),fill=colpal[c(1,2,3,4)])

for(i in c(1,2,3,4)){
  hmp1<-read.table(hmplist[i],fill=T,header=T)
  summy<-apply(log(hmp1[,5:ncol(hmp1)]+0.01),2,mean)
  plot(1:length(summy),summy,type='l',col=colpal[i],lwd=2,ylim=c(-4.7,log(ymax)),ylab='log(+0.01)');par(new=T)
}

