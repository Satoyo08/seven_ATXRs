source('custom_functions.r')
#
# read alpha_screening population and motif DBs
Jan23<-read.table('../data/AlphaScreen/ATX3_VS_1490TFs.csv',sep=",",header=TRUE,quote="\"") # == Sup table 3
length(unique(Jan23$AGI)) # there are some redundant entry within 1490 TFs, Note AT5G06960 (OBF5/TGA5) has two entry, one with ATX3/WGE =10.4, another with 3.9, the latter having a large value for WGE.

# compare alphascreening and ame
alpha_score<-Jan23
motif_enrich<-read.table('../data/AlphaScreen/ame.tsv',sep='\t',header=T) #== Sup table 4

DF<-merge(motif_enrich,alpha_score,by.x='motif_ID',by.y='AT_ID',all=T)

cpy(unique(sort(DF$motif_ID))) # background reference list for GO analysis
cpy(DF)
dev.off();par(mar=c(2,2,0.2,0.2))
plot(-log10(DF$adj_p.value),DF$X.1,xlim=c(0,max(-log10(DF[DF$adj_p.value>0,]$adj_p.value),na.rm=T)),
     ylim=c(0,max(DF$X.1,na.rm=T)),col=ifelse(DF$adj_p.value<0.01,'tomato','black'),pch=ifelse(DF$WGE>quantile(DF$WGE,0.8,na.rm=TRUE),1,16))
cor.test(-log10(DF$adj_p.value),DF$X.1,method='spearman') # rho=-0.077

HIGH_alpha<-na.omit(DF[DF$X.1>6,])
text(-log10(HIGH_alpha$adj_p.value),HIGH_alpha$X.1,HIGH_alpha$NAME,pos=4,col=ifelse(HIGH_alpha$adj_p.value<0.01,'tomato','black'),cex=0.5)
HIGH_p<-DF[-log10(DF$adj_p.value)>250,]
text(-log10(HIGH_p$adj_p.value),HIGH_p$X.1,HIGH_p$NAME,pos=2,col=ifelse(HIGH_p$adj_p.value<0.01,'tomato','black'),cex=0.5)
legend('top',legend=c('adjusted_p-value < 0.01', '>= 0.01','WGE > upper 20%','=< upper 20%'),
       fill=c('tomato','black','white','grey60'),border=c('tomato','black','grey60','grey60'))

dev.off();par(mar=c(2,2,0.2,0.2),mfrow=c(2,1))
hist(DF$X.1,breaks=200,xlim=c(0,max(DF$X.1,na.rm=T)),border='grey30',col='grey30',main='')
hist(-log10(DF$adj_p.value),breaks=200,xlim=c(0,max(-log10(DF[DF$adj_p.value>0,]$adj_p.value),na.rm=T)),ylim=c(0,30),border='grey30',col='grey30')

