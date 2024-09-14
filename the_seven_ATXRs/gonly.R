gonly<-function(X){ # for a given dataframe, match AtID with araport 11 annotation and return a subdataframe containing "protein coding genes"
  arapo<-read.table("/Users/Satoyo/Desktop/temp/論文作業/論文figs/rename_later/Arabidopsis_H3K4me1/data/refs/araport11_all_sorted.bed",header=F) #V6; 1=pc,2=pseudo,3=TE,4=noncoding,5=noveltranscrived
  colnames(arapo)<-c("Chrr","st","en","ID","d","type","V7")
  X1<-merge(X,arapo,by="ID")
  XX<-(X1[!duplicated(X1$ID),])
  X1<-subset(XX,XX$Chrr!="C"&XX$Chrr!="M"&XX$type==1) 
  return(X1)
}
