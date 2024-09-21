setwd("/Users/Satoyo/Desktop/temp/論文作業/論文figs/rename_later/Arabidopsis_H3K4me1")

ATX1_st<-read.csv("data/Figure1/ATX1.csv",header=TRUE);ATX1_len<-1062
ATX2_st<-read.csv("data/Figure1/ATX2.csv",header=TRUE);ATX2_len<-1083
ATX3_st<-read.csv("data/Figure1/ATX3.csv",header=TRUE);ATX3_len<-1018
ATX4_st<-read.csv("data/Figure1/ATX4.csv",header=TRUE);ATX4_len<-1027
ATX5_st<-read.csv("data/Figure1/ATX5.csv",header=TRUE);ATX5_len<-1043
ATXR7_st<-read.csv("data/Figure1/ATXR7.csv",header=TRUE);ATXR7_len<-1423
ATXR3_st<-read.csv("data/Figure1/ATXR3.csv",header=TRUE);ATXR3_len<-2335


setwd("/Users/Satoyo/Desktop/on_enzymes/protein_structure/domain_annotation/")
STRUCTURE<-list(ATX1_st,ATX2_st,ATX3_st,ATX4_st,ATX5_st,ATXR7_st,ATXR3_st);length<-c(ATX1_len,ATX2_len,ATX3_len,ATX4_len,ATX5_len,ATXR7_len,ATXR3_len)
All_st<-rbind(ATX1_st,ATX2_st,ATX3_st,ATX4_st,ATX5_st,ATXR7_st,ATXR3_st)
outnames<-c("ATX1","ATX2","ATX3","ATX4","ATX5","ATXR7","ATXR3")
# set colors - as before
pallets<-as.data.frame(c(brewer.pal(8, "Set1"),brewer.pal(11, "Paired")[c(9,11,1)]))
rownames(pallets)<-unique(All_st$text)[c(7,1,2,3,4,9,5,6,8,10,11)];colnames(pallets)<-c("colors")
pal<-read.table("../color_code",comment.char = "")
pallets<-data.frame(pal$V2)
colnames(pallets)[1]<-'colors'
rownames(pallets)<-pal$V1
#make chimera script for coloring

for(i in 1:length(length)){
  strline<-data.frame("color #",i,"/A:",STRUCTURE[[i]]$start,"-",STRUCTURE[[i]]$end," ",as.character(pallets[as.character(STRUCTURE[[i]]$text),]))
  write.table(strline,file=paste(outnames[i],"_chimera_colors.txt",sep=""),sep="",quote=F,col.names = F,row.names = F)
}

# 色違うじゃん
# 結局 affinity desinger 上で色つけた
