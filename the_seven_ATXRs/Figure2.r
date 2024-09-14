
# This Script train , save and visualize randam forest models
# 'ChIP14' = dataset shown in Fig2. ChIP15'=biological replicates shown in Supplementry Fig
 

# Note that each training results in slightly different models, due to the random numbers that are part of the modeling process.
# Therefore the model you generate from running this script will not match the models provided in /rep1_models and /rep2_models exactly, but be similar overall.

source("RF_functions.r")

# ------- rep1 --------
## read data for rep1 
savedir='..data/RandomForest/rep1_models'

TSS<-read.table('../data/RandomForest/ChIP14_RPM_TSS.txt',header=T)
TES<-read.table('../data/RandomForest/ChIP14_RPM_TES.txt',header=T)
gb<-read.table('../data/RandomForest/ChIP14_RPM.txt',header=T)
A<-merge(TSS,TES,by='ID');B<-gonly(merge(A,gb,by='ID'));AAA<-gonly(merge(A,gb,by='ID'));B[,212:221]<-1000*AAA[,212:221]/(AAA$en-AAA$st)
C<-B[,c(1,grep('TAG',colnames(B)))];tags<-C[,grep('v',colnames(C),invert=T)];nams<-colnames(tags[18:ncol(tags)])
colnames(tags)<-c('ID',paste(nams,'TSS',sep='_'),paste(nams,'TES',sep='_'),paste(nams,'gb_RPKM',sep='_'))



### do RF
# 1. ATX1_TSS     
S<-tags[order(tags$WT.TAG_TSS-tags$ATX1.TAG_TSS),];posID<-S[1:3000,1];
ATX1_TSS<-repeat_train(AL,posID,'remove nothing')
savefilename="ATX1_TSS"
save(ATX1_TSS,file=paste(savedir,savefilename,sep=''))

# 2. ATX2_TSS
S<-tags[order(tags$WT.TAG_TSS-tags$ATX2.TAG_TSS),];posID<-S[1:3000,1];
ATX2_TSS<-repeat_train(AL,posID,'remove nothing')
savefilename="ATX2_TSS"
save(ATX2_TSS,file=paste(savedir,savefilename,sep=''))

# 3. ATX3_TSS
S<-tags[order(tags$WT.TAG_TSS-tags$ATX3.TAG_TSS),];posID<-S[1:3000,1];
ATX3_TSS<-repeat_train(AL,posID,'remove nothing')
savefilename="ATX3_TSS"
save(ATX3_TSS,file=paste(savedir,savefilename,sep=''))

# 4. ATX4_TSS
S<-tags[order(tags$WT.TAG_TSS-tags$ATX4.TAG_TSS),];posID<-S[1:3000,1];
ATX4_TSS<-repeat_train(AL,posID,'remove nothing')
savefilename="ATX4_TSS"
save(ATX4_TSS,file=paste(savedir,savefilename,sep=''))

# 5. ATX5_TSS
S<-tags[order(tags$WT.TAG_TSS-tags$ATX5.TAG_TSS),];posID<-S[1:3000,1];
ATX5_TSS<-repeat_train(AL,posID,'remove nothing')
savefilename="ATX5_TSS"
save(ATX5_TSS,file=paste(savedir,savefilename,sep=''))

# 6. ATXR3_TSS
S<-tags[order(tags$WT.TAG_TSS-tags$ATXR3.TAG_TSS),];posID<-S[1:3000,1];
ATXR3_TSS<-repeat_train(AL,posID,'remove nothing')
savefilename="ATXR3_TSS"
save(ATXR3_TSS,file=paste(savedir,savefilename,sep=''))

# 7. ATXR7_TSS
S<-tags[order(tags$WT.TAG_TES-tags$ATXR7.TAG_TES),];posID<-S[1:3000,1];
ATXR7_TES<-repeat_train(AL,posID,'remove nothing')
savefilename="ATXR7_TES"
save(ATXR7_TES,file=paste(savedir,savefilename,sep=''))

### 4.2 plot RF ---------
model_list<-list.files(savedir,full.names = T)
par(mar=c(0.5,2,0.5,0),mfrow=c(7,2))
for(i in 7:1){
  model<-loadRData(model_list[i])
  plot_RFmodel(model,pallets,sort,col_key,space)
  #mtext(side=3,basename(model_list[i]),line = -1)
}

dev.off();par(mar=c(2,2,0.5,0.5),mfrow=c(1,7))
for(i in 1:7){
  model<-loadRData(model_list[i])
  myROCplot(model,error_FPR_position,ref)
  #mtext(side=3,basename(model_list[i]),line = -1)
}

### are those bound regions open or closed -------
S<-tags[order(tags$WT.TAG_TSS-tags$ATX3.TAG_TSS),];ATX3<-S[1:3000,1];
S<-tags[order(tags$WT.TAG_TSS-tags$ATX4.TAG_TSS),];ATX4<-S[1:3000,1];
S<-tags[order(tags$WT.TAG_TSS-tags$ATX5.TAG_TSS),];ATX5<-S[1:3000,1];

boxplot(AL$TSS_H1[AL$ID %in% ATX3],AL$TSS_H1[AL$ID %in% ATX4],AL$TSS_H1[AL$ID %in% ATX5],AL$TSS_H1,outline=F,notch=T,col=c(pallets[c(7,8,9)],'gray80'))


# ------- rep2 --------
# read data for rep2 
savedir='..data/RandomForest/rep2_models'
TSS<-read.table('../data/RandomForest/ChIP15_TAGS_RPM_TSS.txt',header=T)
TES<-read.table('../data/RandomForest/ChIP15_TAGS_RPM_TES.txt',header=T)
gb<-read.table('../data/RandomForest/ChIP15_TAGS_RPM.txt',header=T)
tags<-gonly(merge(TSS,TES,by='ID'))[1:17]
nams<-colnames(tags[2:9])
colnames(tags)<-c('ID',paste(nams,'TSS',sep='_'),paste(nams,'TES',sep='_'))


# ATX3_TSS
S<-tags[order(tags$WT.rep2.x_TSS-tags$ATX3.TAG.rep2.x_TSS),];posID<-S[1:3000,1];
ATX3_TSS<-repeat_train(AL,posID,'remove nothing')
savefilename="ATX3_TSS"
save(ATX3_TSS,file=paste(savedir,savefilename,sep=''))

# ATX4_TSS
S<-tags[order(tags$WT.rep2.x_TSS-tags$ATX4.TAG.rep2.x_TSS),];posID<-S[1:3000,1];
ATX4_TSS<-repeat_train(AL,posID,'remove nothing')
savefilename="ATX4_TSS"
save(ATX4_TSS,file=paste(savedir,savefilename,sep=''))

# ATX5_TSS
S<-tags[order(tags$WT.rep2.x_TSS-tags$ATX5.TAG.rep2.x_TSS),];posID<-S[1:3000,1];
ATX5_TSS<-repeat_train(AL,posID,'remove nothing')
savefilename="ATX5_TSS"
save(ATX5_TSS,file=paste(savedir,savefilename,sep=''))

# ATXR3_TSS
S<-tags[order(tags$WT.rep2.x_TSS-tags$ATXR3.TAG.rep2.x_TSS),];posID<-S[1:3000,1];
ATXR3_TSS<-repeat_train(AL,posID,'remove nothing')
savefilename="ATXR3_TSS"
save(ATXR3_TSS,file=paste(savedir,savefilename,sep=''))
