# This Script reads the output files from the trained SVM models and visualizes them
# The script for training can be found at ../scripts/SVM_training.ipyenv

### 5.1 output ROC
library(Biostrings);library(qgraph)
source('custom_functions.r')
# for rep1
ROC1<-plot_ROC_fromTPRFPR('../data/SVM/rep1/SVM_outputs/ATX1_TSS_6mer_FPRTPR') 
ROC2<-plot_ROC_fromTPRFPR('../data/SVM/rep1/SVM_outputs/ATX2_TSS_6mer_FPRTPR') 
ROC3<-plot_ROC_fromTPRFPR('../data/SVM/rep1/SVM_outputs/ATX3_TSS_6mer_FPRTPR')
ROC4<-plot_ROC_fromTPRFPR('../data/SVM/rep1/SVM_outputs/ATX4_TSS_6mer_FPRTPR') 
ROC5<-plot_ROC_fromTPRFPR('../data/SVM/rep1/SVM_outputs/ATX5_TSS_6mer_FPRTPR') 
ROCR3<-plot_ROC_fromTPRFPR('../data/SVM/rep1/SVM_outputs/ATXR3_TSS_6mer_FPRTPR') 
ROCR7<-plot_ROC_fromTPRFPR('../data/SVM/rep1/SVM_outputs/ATXR7_TES_6mer_FPRTPR') 
# for rep2
ROC3<-plot_ROC_fromTPRFPR('../data/SVM/rep2/SVM_outputs/ChIP15_ATX3_FPRTPR_Mar31') 
ROC4<-plot_ROC_fromTPRFPR('../data/SVM/rep2/SVM_outputs/ChIP15_ATX4_FPRTPR_Mar31') 
ROC5<-plot_ROC_fromTPRFPR('../data/SVM/rep2/SVM_outputs/ChIP15_ATX5_FPRTPR_Mar31') 
ROCR3<-plot_ROC_fromTPRFPR('../data/SVM/rep2/SVM_outputs/ChIP15_ATXR3_FPRTPR_Mar31') 


ROC1$auc_score # get ROC score. for ROC1, for example

### 5.2 load weights --------
# rep 1
weightfiles<-list.files('../data/SVM/rep1/SVM_outputs',pattern='weigths',full.name=TRUE);weightfiles
ST1<-read.table(weightfiles[1],header=F)
for(i in 1:7){
  weights1<-read.table(weightfiles[i],header=F);weights1<-add_mean(weights1);
  ST1[,1+i]<-weights1$V7
  colnames(ST1)[1+i]<-substr(basename(weightfiles[i]),1,nchar(basename(weightfiles[i]))-17)
}

# rep 2 
weightfiles<-list.files('../data/SVM/rep2/SVM_outputs',pattern='weights_Mar31',full.name=TRUE);weightfiles
ST2<-read.table(weightfiles[1],header=F)[,1:5]
for(i in 1:4){
  weights1<-read.table(weightfiles[i],header=F);weights1<-add_mean(weights1);
  ST2[,1+i]<-weights1$V7
  colnames(ST2)[1+i]<-substr(basename(weightfiles[i]),1,nchar(basename(weightfiles[i]))-14)
}

#Oya 2022
weightfiles<-list.files('/Users/Satoyo/Desktop/temp/論文作業/論文figs/rename_later/Arabidopsis_H3K4me1/data/Figure4/weights',full.name=TRUE);weightfiles # replace the path with https://github.com/Satoyo08/Arabidopsis_H3K4me1/tree/main/data/Figure4/weights
ST3<-read.table(weightfiles[1],header=F)[,1:5]
for(i in 1:4){
  weights1<-read.table(weightfiles[i],header=F);weights1<-add_mean(weights1);
  ST3[,1+i]<-weights1$V7
  colnames(ST3)[1+i]<-substr(basename(weightfiles[i]),1,nchar(basename(weightfiles[i]))-8)
}
head(ST3)
# for supplementary Fig, compare weights
par(mfrow=c(2,4),mar=c(2,2,0.2,0.2))
m1<-ST1$ATX1;m2<-ST3$ChIP7_ATX1;plot(m1,m2,pch=16,cex=0.6);text(min(m1),max(m2)*0.9,paste("r = ",round(cor.test(m1,m2,method="pearson")$estimate,digits=2)),pos=4,cex=1.5)
m1<-ST1$ATX2;m2<-ST3$ChIP7_ATX2;plot(m1,m2,pch=16,cex=0.6);text(min(m1),max(m2)*0.9,paste("r = ",round(cor.test(m1,m2,method="pearson")$estimate,digits=2)),pos=4,cex=1.5)
plot(1,1,type='n');plot(1,1,type='n',xaxt='n')
m1<-ST1$ATX3;m2<-ST2$ChIP15_ATX3;plot(m1,m2,pch=16,cex=0.6);text(min(m1),max(m2)*0.9,paste("r = ",round(cor.test(m1,m2,method="pearson")$estimate,digits=2)),pos=4,cex=1.5)
m1<-ST1$ATX4;m2<-ST2$ChIP15_ATX4;plot(m1,m2,pch=16,cex=0.6);text(min(m1),max(m2)*0.9,paste("r = ",round(cor.test(m1,m2,method="pearson")$estimate,digits=2)),pos=4,cex=1.5)

m1<-ST1$ATX5;m2<-ST2$ChIP15_ATX5;plot(m1,m2,pch=16,cex=0.6,xlim=c(-0.006,0.009),ylim=c(-0.007,0.012));text(min(m1),max(m2)*0.9,paste("r = ",round(cor.test(m1,m2,method="pearson")$estimate,digits=2)),pos=4,cex=1.5)
m1<-ST1$ATXR3;m2<-ST2$ChIP15_ATXR3;plot(m1,m2,pch=16,cex=0.6);text(min(m1),max(m2)*0.9,paste("r = ",round(cor.test(m1,m2,method="pearson")$estimate,digits=2)),pos=4,cex=1.5)

### 5.3 plot string --------

replace_features<-function(feature_list){　# with alpha_screening features
  border_colors=rep('gray80',length(feature_list)) 
  border_colors<-replace(border_colors,grep('AAGAGA|GAGAGA|AAGAGG|GAGAGG|AGAGAA|GGAGAA|AGAGAG|GGAGAG',feature_list),'gray30');border_colors  # RAGAGR or RGAGAR
  border_colors<-replace(border_colors,grep('TCTCTT|TCTCTC|CCTCTT|CCTCTC|TTCTCT|TTCTCC|CTCTCT|CTCTCC',feature_list),'gray30');border_colors  # RAGAGR or RGAGAR
  border_colors<-replace(border_colors,grep('AAACCC|AACCCT|ACCCAA|ACCCTA',feature_list),"#800080") ;border_colors  #telobox
  border_colors<-replace(border_colors,grep('AAGCCC|AGGCCC|GGCCCA|AGCCCA|GCCCAA|GCCCAT|GGGCTT|GGGCCT|TGGGCC|TGGGCT|TTGGGC|ATGGGC|AATGGG|CCCATT',feature_list),pallets[3]);border_colors #my own consensus for ARGCCCAWT 
  border_colors<-replace(border_colors,grep('GCGGCG|CGCCGC|CCGGCG|CGCCGG|CGGCGG|CCGCCG|CGGCGA|TCGCCG|TCGGCG|CGCCGA',feature_list),pallets[1]);border_colors  # from ATX3 string, SCGGCGR. also RAP2.11 site
  border_colors<-replace(border_colors,grep('GTCGTC|GACGAC|TCGTCG|CGACGA|CGTCGT|ACGACG',feature_list),pallets[5]);border_colors  # from ATX3 string, TCGTCGTC
  border_colors<-replace(border_colors,grep('TCGGAA|TTCCGA|ATCGGA|TCCGAT|AATCGG|CCGATT|GAATCG|CGATTC|AATCGA|TCGATT',feature_list),pallets[8]);border_colors  # from ATX3 string, TCYGATTC
  border_colors<-replace(border_colors,grep('TCTGCA|CTCTGC|TCTCTG|CTGCAA|TGCAGA',feature_list),pallets[10]);border_colors  # from ATX2 string, JMJ28
  return(border_colors)
}

# plot legends---
dev.off()
plot(rep(1,7),1:7,pch=16,cex=3,col=c('gray30','#800080',pallets[c(3,1,5,8,10)]))
text(rep(1,7),1:7,pos=4,c('RAGAR','telobox','ARGCCCAWT','SCGGCGR=RAP2.11','TCGTCGTC','TCYGATTC','JMJ28'))

# plot qgraphs--
dev.off()
height=64/25.4;width=112/25.4;pointsize=12
library(Biostrings);library(qgraph)

# for_rep1
file_prefix_list=c('ATX1_rep1','ATX2_rep1','ATX3_rep1','ATX4_rep1',"ATX5_rep1","ATXR3_rep1","ATXR7_rep1")
# for_rep2
file_prefix_list=c('ATX3_ChIP15_Mar31','ATX4_ChIP15_Mar31','ATX5_ChIP15_Mar31','ATXR3_ChIP15_Mar31')



for(i in 1:length(weightfiles)){
  file_prefix=file_prefix_list[i]  
  W<-add_mean(read.table(weightfiles[i]));weightfiles[i]
  Sorted_W<-W[order(-W[,7]),]
  for(j in 1:2){
    if(j==1){
      # positive k-mer
      cairo_pdf(file=paste(file_prefix,"pos_cairo.pdf",sep='_'),height=height,width=width,pointsize=pointsize)
      mer<-data.frame(Sorted_W[1:60,c(1,7)],rownames(Sorted_W[1:60,c(1,7)]));colnames(mer)<-c('kmer','weights','index')
    }else{
      cairo_pdf(file=paste(file_prefix,"neg_cairo.pdf",sep='_'),height=height,width=width,pointsize=pointsize)
      mer<-data.frame(Sorted_W[nrow(W):(nrow(W)-59),c(1,7)],rownames(Sorted_W[nrow(W):(nrow(W)-59),c(1,7)]));colnames(mer)<-c('kmer','weights','index')
    }
    feature_list<-as.character(mer$kmer);feature_list_s<-sort(feature_list)
    DNAS<-list();for(i in 1:length(feature_list)){DNAS[i]<-DNAString(feature_list_s[i])}
    RC_DNAS<-list();for(i in 1:length(feature_list)){RC_DNAS[i]<-reverseComplement(DNAS[[i]])}
    M<-kmer_list_to_matrix_dir(feature_list);dist_mi<-M[[1]];X<-M[[2]] 
    
    border_colors<-replace_features(feature_list)
    charsize<-2+abs(mer[order(feature_list),2]*180);charsize
    #draw qgraph!
    Q<-qgraph(dist_mi, layout='spring',vsize=ifelse(charsize <= 3.5, 3.5, charsize),
              color='gray80',border.color=border_colors[order(feature_list)],
              border.width=2.5,labels=feature_list_s,label.cex=2.5,edge.color='forestgreen')
    mtext(side=3,basename(weightfiles[i]),line=-1.5)
    dev.off()
  }
}

### 5.4 weights Heatmaps
rankST<-ST1 
for(i in 1:7){rankST[,i+1]<-rank(ST1[,i+1])} # negative weights = rank high (1~)
rankTF<-ST1
for(i in 2:8){rankTF[,i]<-c(rankST[,i]<61 | rankST[,i]>(4096-60))}
for(i in 2:8){rankTF[,i]<-c(rankST[,i]>(4096-60))}

ST1[,9]<-apply(rankTF[,2:8],1,sum)
feature_list<-ST1[ST1$V9>0,1]
border_colors<-replace_features(feature_list)
hmp_autoclust<-heatmap(as.matrix(ST1[ST1$V9>0,2:8]),Colv=NA,scale='col',RowSideColors=border_colors,labRow=F)#ST1[ST1$V9>0,1]) 

### 5.5 motif distribution over TSS ------
setwd("../data/ATX3_motif_search_in_TSS")
# read the bedfile, merge with original 3000 gene list, sequentially numbered. Use the sequential number and location to make polygon images
system('grep ">" ATX3_TSS_bound_TSS_150_300.fasta | sed -e "s/>//" > posheaders')
system('grep ">" random_3000_ATX3_unbound.fasta | sed -e "s/>//" > negheaders')
all_posid<-read.table('posheaders')
all_negid<-read.table('negheaders')
all_posid$sequancial<-1:nrow(all_posid);nrow(all_posid)
all_negid$sequancial<-1:nrow(all_negid);nrow(all_negid)

# load bedfiles representing the motif location- these files were generated using seqkit;

# posfa=ATX3_TSS_bound_TSS_150_300.fasta
# negfa=ATX3_TSS_unbound_TSS_150_300.fasta
# seqkit sample -n 3077 $negfa > random_3000_ATX3_unbound.fasta # -n and output line number does not exactly match. -n 3077 this achieves 3001 sequence
# negsample=random_3000_ATX3_unbound.fasta
# seqkit locate -p 'TCGTCGTC' -m 0 $posfa > ATX3_TCGTCGTC.bed
# seqkit locate -p 'TCGTCGTC' -m 0 $negsample > ATX3_unbound_TCGTCGTC.bed
# 
# seqkit locate -p 'TCCGATTC' -p 'TCTGATTC' -m 0 $posfa > ATX3_TCYGATTC.bed
# sequin locate -p 'TCCGATTC' -p 'TCTGATTC' -m 0 $negsample > ATX3_unbound_TCYGATTC.bed
# 
# # SCGGCGR (S = G or C, R = A or G),
# seqkit locate -p 'GCGGCGA' -p 'GCGGCGG' -p 'CCGGCGA' -p 'CCGGCGG' -m 0 $posfa > ATX3_SCGGCGR.bed
# seqkit locate -p 'GCGGCGA' -p 'GCGGCGG' -p 'CCGGCGA' -p 'CCGGCGG' -m 0 $negsample > ATX3_unbound_SCGGCGR.bed
# 
# # RGCCCAW (R = A or G, W = A or T) 
# seqkit locate -p 'AGCCCAA' -p 'AGCCCAT' -p 'GGCCCAA' -p 'GGCCCAT' -m 0 $posfa > ATX3_RGCCCAW.bed
# seqkit locate -p 'AGCCCAA' -p 'AGCCCAT' -p 'GGCCCAA' -p 'GGCCCAT' -m 0 $negsample > ATX3_unbound_RGCCCAW.bed

bound_bed<-read.table("ATX3_RGCCCAW.bed",header=T)
unbound_bed<-read.table("ATX3_unbound_RGCCCAW.bed",header=T)

POSmergedDF<-merge(all_posid,bound_bed,by.x='V1',by.y='seqID',all.x=TRUE)
NEGmergedDF<-merge(all_negid,unbound_bed,by.x='V1',by.y='seqID',all.x=TRUE)


PhistDF<-hist(POSmergedDF$start,xlim=c(0,451),breaks=10)
NhistDF<-hist(NEGmergedDF$start,xlim=c(0,451),breaks=10)

#colors; TGCTGCTC #f5ca51 TCYGATTC #ABDDA4, SCGGCGR #bd0954, RGCCCAW  #FDAD62
# EPS size 245 x 252

clr="#FDAD62"
dev.off();par(mfrow=c(3,1));par(mar=c(0.2,2,0.2,0.2))
plot(1:length(NhistDF$counts),PhistDF$counts,col=clr,lwd=2,type='l',xaxt='n',ylim=c(0,max(PhistDF$counts)));par(new=T)
plot(1:length(NhistDF$counts),NhistDF$counts,col="grey60",type='l',lwd=2,xaxt='n',yaxt='n',ylim=c(0,max(PhistDF$counts)))

plot(1,1,type='n',xlim=c(0,451),ylim=c(0,3000),xaxt='n',yaxt='n',xlab='',ylab='')
abline(v=150,lty=2)
rect(POSmergedDF$start,POSmergedDF$sequancial-1,POSmergedDF$end,POSmergedDF$sequancial,col=clr,border=clr)

plot(1,1,type='n',xlim=c(0,451),ylim=c(0,3000),xaxt='n',yaxt='n',xlab='',ylab='')
abline(v=150,lty=2)
rect(NEGmergedDF$start,NEGmergedDF$sequancial-1,NEGmergedDF$end,NEGmergedDF$sequancial,col=clr,border=clr)

