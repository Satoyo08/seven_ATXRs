
refg_g<-read.table("../araport_genes_and_polyA_gonly.bed",header=F,fill=TRUE);colnames(refg_g)[10]<-'ID';head(refg_g)
# set TTS of the gene
refg_g$genes_3p<-NA
refg_g[refg_g$V6=='+',]$genes_3p<-refg_g[refg_g$V6=='+',]$V9
refg_g[refg_g$V6=='-',]$genes_3p<-refg_g[refg_g$V6=='-',]$V8

# set PAC start site 
refg_g$PAC_5<-NA
refg_g[refg_g$V6=='+',]$PAC_5<-refg_g[refg_g$V6=='+',]$V2
refg_g[refg_g$V6=='-',]$PAC_5<-refg_g[refg_g$V6=='-',]$V3

# the distance
refg_g$dist<-refg_g$genes_3p-refg_g$PAC_5
# numbering of polyA, form 3' to 5' of gene
refg_g$rank<-unlist(tapply(abs(refg_g$dist),droplevels(as.factor(refg_g$V11)),rank));head(refg_g)


# dataset1 ------------------
filelist<-list.files(path='../polyA_dataset1/',pattern='bed',full.names = TRUE);filelist
RNA8list<-list()
for(i in 1:length(filelist)){
  R<-read.table(filelist[i],header=F)
  C<-R[-grep('t',R$V1),]
  name<-substr(basename(filelist[i]),1,nchar(basename(filelist[i]))-25)
  DF<-merge(C,refg_g,by.x="V10",by.y="V4",all.x=TRUE)
  print(c(nrow(DF),nrow(C)))
  RNA8list[[i]]<-DF
  names(RNA8list)[i]<-name
}

atxr7<-rbind(RNA8list$atxr7_1,RNA8list$atxr7_2,RNA8list$atxr7_3);nrow(atxr7)
wt<-rbind(RNA8list$wt_1,RNA8list$wt_2,RNA8list$wt_3);nrow(wt)
atxr7tab<-table(atxr7$rank)
wttab<-table(wt$rank)

mean(atxr7$rank,na.rm=T);mean(wt$rank,na.rm=T) #average position of mapped polyA site as a rank from 3' most polyA site


# datasets2 ------------------
filelist<-list.files(path='../polyA_dataset2/',full.names = TRUE,pattern='109');filelist
vernRNA<-data.frame()
for(i in 1:length(filelist)){
  R<-read.table(filelist[i],header=F)
  C<-R[-grep('t',R$V1),]
  name<-substr(basename(filelist[i]),15,nchar(basename(filelist[i]))-28)
  DF<-merge(C,refg_g,by.x="V10",by.y="V4",all.x=TRUE)
  print(c(name,nrow(DF),nrow(C)))
  vernRNA<-rbind(vernRNA,DF)
}
vernATXR7s<-vernRNA

filelist<-list.files(path='../polyA_dataset2/',full.names = TRUE,pattern='122');filelist
vernRNA<-data.frame()
for(i in 1:length(filelist)){
  R<-read.table(filelist[i],header=F)
  C<-R[-grep('t',R$V1),]
  name<-substr(basename(filelist[i]),15,nchar(basename(filelist[i]))-28)
  DF<-merge(C,refg_g,by.x="V10",by.y="V4",all.x=TRUE)
  print(c(name,nrow(DF),nrow(C)))
  vernRNA<-rbind(vernRNA,DF)
}
vernWTs<-vernRNA

vern_atxr7_tab<-table(vernATXR7s$rank)
vern_WT_tab<-table(vernWTs$rank)

mean(vernATXR7s$rank,na.rm=T);mean(vernWTs$rank,na.rm=T)
