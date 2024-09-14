

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

### 5.2 compare weights --------
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

# -> ST1 and ST2 were combined and provided as Supplementary table with the manuscript
# STable<-merge(ST1,ST2,by="V1")
# colnames(STable)[c(1,9:12)]<-c("6-mer","rep2_ATX3","rep2_ATX4","rep2_ATX5","rep2_ATXR3")
# write.table(STable,file='../data/SVM/Combined_SVM_weights.txt',sep='\t',quote=F,row.names = F)

### 5.3 plot string --------

#dev.off();par(mar=c(0,0,0,0),mfrow=c(1,2))
#for(i in 1:7){
W<-add_mean(read.table(weightfiles[i]))
Sorted_W<-W[order(-W[,7]),]

replace_features_old<-function(feature_list){
  border_colors=rep('gray80',length(feature_list)) 
  pal<-c('gray30',pallets[7],pallets[4],'purple')
  border_colors<-replace(border_colors,grep('AAGAGA|GAGAGA|AAGAGG|GAGAGG|AGAGAA|GGAGAA|AGAGAG|GGAGAG',feature_list),pal[1]);border_colors  # RAGAGR or RGAGAR
  border_colors<-replace(border_colors,grep('TCTCTT|TCTCTC|CCTCTT|CCTCTC|TTCTCT|TTCTCC|CTCTCT|CTCTCC',feature_list),'gray60');border_colors  # RAGAGR or RGAGAR
  border_colors<-replace(border_colors,grep('AAACCC|AACCCT|ACCCAA|ACCCTA',feature_list),pal[2]) ;border_colors  #telobox
  border_colors<-replace(border_colors,grep('AAGCCC|AGGCCC|GGCCCA|AGCCCA|GCCCAA|GCCCAT|GGGCTT|GGGCCT|TGGGCC|TGGGCT|TTGGGC|ATGGGC|AATGGG|CCCATT',feature_list),pal[3]);border_colors #my own consensus for ARGCCCAWT 
  border_colors<-replace(border_colors,grep('TATAAA|TATATA|ATATAA|ATAAAT|TAAATA|ATATAT|TTATAA|TTATAT',feature_list),pal[4]);border_colors  # TATA box from yamamoto 
  border_colors<-replace(border_colors,grep('GCGGCG|CGCCGC|CCGGCG|CGCCGG|CGGCGG|CCGCCG|CGGCGA|TCGCCG|TCGGCG|CGCCGA',feature_list),pallets[2]);border_colors  # from ATX3 string
  border_colors<-replace(border_colors,grep('GTCGTC|GACGAC|TCGTCG|CGACGA|CGTCGT|ACGACG',feature_list),pallets[4]);border_colors  # from ATX3 string, TCGTCGTC
  border_colors<-replace(border_colors,grep('GTCGTC|GACGAC|TCGTCG|CGACGA|CGTCGT|ACGACG',feature_list),pallets[6]);border_colors  # from ATX3 string, TCGTCGTC
  border_colors<-replace(border_colors,grep('TCGGAA|TTCCGA|ATCGGA|TCCGAT|AATCGG|CCGATT|GAATCG|CGATTC|AATCGA|TCGATT',feature_list),pallets[8]);border_colors  # from ATX3 string, TCGTCGTC
  border_colors<-replace(border_colors,grep('TCTGCA|CTCTGC|TCTCTG|CTGCAA|TGCAGA',feature_list),pallets[10]);border_colors  # from ATX2 string
  return(border_colors)
}

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
  #alpha_screening features
  # border_colors<-replace(border_colors,grep('CCCTTA|CCCTAA|TACCCTT|TGCCCTT|TAAGGG|TTAGGG',feature_list),pallets[1]);border_colors  #ASR3
  # border_colors<-replace(border_colors,grep('GCCGGC|GCCGTC',feature_list),pallets[3]);border_colors  #RAP2
  # border_colors<-replace(border_colors,grep('TGACGT|GACGTC|ACGTCA',feature_list),pallets[5]);border_colors  #TGA6/OBF/ATTGA2
  # border_colors<-replace(border_colors,grep('CACGTG',feature_list),pallets[11]);border_colors  #bLHL
  # border_colors<-replace(border_colors,grep('TTGACT|TTGACC|AGTCAA|GGTCAA',feature_list),'#800080');border_colors  #WRKY40
  return(border_colors)
}


plot(rep(1,7),1:7,pch=16,cex=3,col=c('gray30','#800080',pallets[c(3,1,5,8,10)]),axes=F,xlab='',ylab='')
text(rep(1.1,7),1:7,pos=4,c('RAGAR','telobox','ARGCCCAWT','SCGGCGR=RAP2.11','TCGTCGTC','TCYGATTC','JMJ28'))



#qgraphs--
dev.off()
height=25/25.4;width=45/25.4;pointsize=6
install.packages('extrafont')
library(extrafont)
font_import()  # This imports and registers the system fonts, might take a few minutes
fonts()  # Lists all available fonts... does not work.


file_prefix_list=c('ATX1_ChIP14','ATX2_ChIP14','ATX3_ChIP14','ATX4_ChIP14','ATX5_ChIP14','ATXR3_ChIP14','ATXR7_TES_ChIP14')
for(i in 1:7){
  file_prefix=file_prefix_list[i]  
  W<-add_mean(read.table(weightfiles[i]));weightfiles[i]
  Sorted_W<-W[order(-W[,7]),]
  for(j in 1:2){
    if(j==1){
      # positive k-mer
      pdf(file=paste(file_prefix,"pos2.pdf",sep='_'),height=height,width=width,pointsize=pointsize,fonts='serif')
      mer<-data.frame(Sorted_W[1:60,c(1,7)],rownames(Sorted_W[1:60,c(1,7)]));colnames(mer)<-c('kmer','weights','index')
    }else{
      pdf(file=paste(file_prefix,"neg2.pdf",sep='_'),height=height,width=width,pointsize=pointsize,fonts='serif')
      mer<-data.frame(Sorted_W[nrow(W):(nrow(W)-59),c(1,7)],rownames(Sorted_W[nrow(W):(nrow(W)-59),c(1,7)]));colnames(mer)<-c('kmer','weights','index')
    }
    feature_list<-as.character(mer$kmer);feature_list_s<-sort(feature_list)
    DNAS<-list();for(i in 1:length(feature_list)){DNAS[i]<-DNAString(feature_list_s[i])}
    RC_DNAS<-list();for(i in 1:length(feature_list)){RC_DNAS[i]<-reverseComplement(DNAS[[i]])}
    M<-kmer_list_to_matrix_dir(feature_list);dist_mi<-M[[1]];X<-M[[2]] 
    
    border_colors<-replace_features(feature_list)
    #draw qgraph!
    Q<-qgraph(dist_mi, layout='spring',vsize=1+abs(mer[order(feature_list),2]*300),
              color='gray80',border.color=border_colors[order(feature_list)],
              border.width=2.5,labels=feature_list_s,label.cex=2.5,edge.color='forestgreen')
    mtext(side=3,basename(weightfiles[i]),line=-1.5)
    dev.off()
  }
}

# pdf がaffinity 上でラスター的に見える。R出力のpdfをロードすると出るエラー missing font ZapfDingbats のせいだろうか？
# ppt だとベクター画像



#5.5 weights heatmap
rankST<-ST1 
for(i in 1:7){rankST[,i+1]<-rank(ST1[,i+1])} # negative weights = rank high (1~)
rankTF<-ST1
for(i in 2:8){rankTF[,i]<-c(rankST[,i]<61 | rankST[,i]>(4096-60))}
for(i in 2:8){rankTF[,i]<-c(rankST[,i]>(4096-60))}

ST1[,9]<-apply(rankTF[,2:8],1,sum)
border_colors<-replace_features(feature_list)
hmp_autoclust<-heatmap(as.matrix(ST1[ST1$V9>0,2:8]),Colv=NA,scale='col',RowSideColors=border_colors,labRow=F)#ST1[ST1$V9>0,1]) 

### 5.* 坂本さんに送るATX345 motifs -----
write.table(ST1[,c(1,4,5,6)],file='ATX345_SVM_weights.txt',sep='\t',quote=F,row.names = F)

### select 'replace features' . What TFs are used in the... ----

# read rep2 dataset
weightfiles<-list.files('/Users/Satoyo/data/ChIP15/rep_ATXRs-all/SVM',pattern='weights_Mar31',full.name=TRUE);weightfiles
ST2<-read.table(weightfiles[1],header=F)[,1:4]
for(i in 1:4){
  weights1<-read.table(weightfiles[i],header=F);weights1<-add_mean(weights1);
  ST2[,1+i]<-weights1$V7
  colnames(ST2)[1+i]<-substr(basename(weightfiles[i]),8,nchar(basename(weightfiles[i]))-8)
}

#?ST<-merge(ST1,ST2,by='V1')
ST<-ST2
rankST<-ST
for(i in 1:4){rankST[,i+1]<-rank(ST[,i+1])} # negative weights = rank high (1~)
rankTF<-ST
for(i in 2:5){rankTF[,i]<-c(rankST[,i]<61 | rankST[,i]>(4096-60))}
#for(i in 2:8){rankTF[,i]<-c(rankST[,i]>(4096-60))}
ST[,ncol(ST)+1]<-apply(rankTF[,2:5],1,sum)
features_of_interest<-ST[ST$V12>0,]$V1 # rep2 ATX345とrep1のATX1-R3に共通してhigh weight なfeatures

cpy(table(replace_features(features_of_interest)))
# #3288BD #5E4FA2 #800080 #ABDDA4 #D53E4F #E6F598 #F46D43 #FDAE61 #FFFFBF  gray30  gray80 
# --> comparison revealed ASR3, TGA6/OBF/ATTGA2 are not highly weighted 

rankTF[rankTF$V1 %in% c("TTGACT","TTGACC","AGTCAA","GGTCAA"),] # WRKY40 binding site--- was enriched for ATXR3 binding site... 
## TSS motif vs ATX 分布 HMP ----- 
#fig b のHMPの代わりに、TSS中のmotifの有無と、横にATX たちの分布量とか書いてみたいなあ、
#RGCCCAWに加えてJMJ28があるグループではATX2が多いね、みたいなのが可視化されるといい

kmers<-read.table(file = '/Users/Satoyo/data/ChIP15/rep_ATXRs-all/SVM/frequency_of_kmers_par_genes.csv',sep=',',header=TRUE)

selected_kmers<-c("AAGAGA","GAGAGA","AAGAGG","GAGAGG","AGAGAA","GGAGAA","AGAGAG","GGAGAG","TCTCTT","TCTCTC","CCTCTT","CCTCTC","TTCTCT","TTCTCC","CTCTCT","CTCTCC","AAACCC","AACCCT","ACCCAA","ACCCTA","AAGCCC","AGGCCC","GGCCCA","AGCCCA","GCCCAA","GCCCAT","GGGCTT","GGGCCT","TGGGCC","TGGGCT","TTGGGC","ATGGGC","AATGGG","CCCATT","GCGGCG","CGCCGC","CCGGCG","CGCCGG","CGGCGG","CCGCCG","CGGCGA","TCGCCG","TCGGCG","CGCCGA","GTCGTC","GACGAC","TCGTCG","CGACGA","CGTCGT","ACGACG","TCGGAA","TTCCGA","ATCGGA","TCCGAT","AATCGG","CCGATT","GAATCG","CGATTC","AATCGA","TCGATT","TCTGCA","CTCTGC","TCTCTG","CTGCAA","TGCAGA")
selected_kmers_df<-kmers[,colnames(kmers) %in% selected_kmers]
selected_kmers_df$ID<-substr(kmers$X,1,9)

head(selected_kmers_df)
tag_norm<-data.frame(tags$ID,tags[,2:8]-tags[,2],tags[,17]-tags[,10]);colnames(tag_norm)[c(1,9)]<-c('ID','ATXR7_TES')

mer_tag<-merge(tag_norm,selected_kmers_df,by='ID')
head(mer_tag)
# ATX1-R3 bound と、random ubound. どんな並び順が良いだろう？


S<-mer_tag[order(mer_tag$ATX2.TAG_TSS),]
heatmap(as.matrix(S[,c(4,10:ncol(S))]),Rowv=NA,Colv=NA,scale='row')

