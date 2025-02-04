

# This Script train, save and visualize randam forest models (Fig.2 and Supplementary Fig.3 )
# 'ChIP14' = dataset shown in Fig2. ChIP15'=biological replicates shown in Supplementry Fig.3
 

# Note that each training results in slightly different models due to the randomness in a part of the modeling process.
# Therefore the model you generate from running this script will not exactly match the models provided in /rep1_models and /rep2_models, but should be similar overall.

source("RF_functions.r")

# ------- rep1 --------
## read data for rep1 
savedir='../data/RandomForest/rep1_models/'

TSS<-read.table('../data/RandomForest/ChIP14_RPM_TSS.txt',header=T)
TES<-read.table('../data/RandomForest/ChIP14_RPM_TES.txt',header=T)
gb<-read.table('../data/RandomForest/ChIP14_RPM.txt',header=T)
A<-merge(TSS,TES,by='ID');B<-gonly(merge(A,gb,by='ID'));AAA<-gonly(merge(A,gb,by='ID'));B[,212:221]<-1000*AAA[,212:221]/(AAA$en-AAA$st)
C<-B[,c(1,grep('TAG',colnames(B)))];tags<-C[,grep('v',colnames(C),invert=T)];nams<-colnames(tags[18:ncol(tags)])
colnames(tags)<-c('ID',paste(nams,'TSS',sep='_'),paste(nams,'TES',sep='_'),paste(nams,'gb_RPKM',sep='_'))



### train RF
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

### 4.2 Visualize RF ---------
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
savedir='../data/RandomForest/rep2_models'
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

# data vis
model_list<-list.files(savedir,full.names = T)
par(mar=c(0.5,2,0.5,0),mfrow=c(4,2))
for(i in 4:1){
  model<-loadRData(model_list[i])
  plot_RFmodel(model,pallets,sort,col_key,space)
  #mtext(side=3,basename(model_list[i]),line = -1)
}

dev.off();par(mar=c(2,2,0.5,0.5),mfrow=c(1,4))
for(i in 1:4){
  model<-loadRData(model_list[i])
  myROCplot(model,error_FPR_position,ref)
  #mtext(side=3,basename(model_list[i]),line = -1)
}


# ----- Compare Importance between replicats (Supplimentary Figure 3) ------
modeldir_rep1='../data/RandomForest/rep1_models'
modeldir_rep2='../data/RandomForest/rep2_models'
modeldir_Oya2022='/Users/Satoyo/Desktop/temp/論文作業/論文figs/rename_later/Arabidopsis_H3K4me1/data/Figure3'# Replace the path to https://github.com/Satoyo08/Arabidopsis_H3K4me1/tree/main/data/Figure3

model_list_rep1<-list.files(modeldir_rep1,full.names = T)[c(3:6)]
model_list_rep2<-list.files(modeldir_rep2,full.names = T)
model_list_rep1_ATX12R37<-list.files(modeldir_rep1,full.names = T)[c(1,2,6,7)]
model_list_Oya2022<-list.files(modeldir_Oya2022,full.names = T)[c(1,3,5,6)]

sort<-c(1,2,24,3,66,45,25,4,67,46,26,5,68,47,36,15,78,57,34,13,76,55,33,12,75,54,35,14,77,56,29,8,71,50,31,10,73,52,30,9,72,51,27,6,69,48,42,21,84,63,43,22,85,64,44,23,86,65,37,16,79,58,28,7,70,49,40,19,82,61,41,20,83,62,32,11,74,53,87,93,90,88,94,91,89,95,92,38,17,80,59,39,18,81,60)
col_key<-c(10,9,rep(c(3:6),19),rep(c(3,4,6),3),rep(c(3:6),2))

radial_shift <- function(x, y, angle, distance) {
  angle_rad <- angle * (pi / 180)
  new_x <- x + distance * sin(angle_rad)  
  new_y <- y + distance * cos(angle_rad)  
  return(c(new_x,new_y))
}

dev.off()
RF_correlation_plot_with_label<-function(Rep1_model_path,Rep2_model_path,threshold,jitter_dist,angle_from,angle_to){
  model1<-loadRData(Rep1_model_path)
  model2<-loadRData(Rep2_model_path)
  IMPs1<-apply(model1$Importance,1,mean)[sort]
  IMPs2<-apply(model2$Importance,1,mean)[sort]
  label_indices <- which(IMPs1 > threshold | IMPs2 > threshold)
  sorted_label_indices<-label_indices[order(IMPs2[label_indices])]
  replace_map <- c("CMA601" = "RNAP total CTD",
                   "CMA602" = "RNAP2 PhosS2",
                   "CMA603" = "RNAP2 PhosS5")
  
  
  plot(IMPs1,IMPs2,col=palette[col_key],pch=16,xlab='',ylab='',xlim=c(0,max(c(IMPs1,IMPs2))*1.1),ylim=c(0,max(c(IMPs1,IMPs2))*1.1))
  dist<-jitter_dist+runif(length(label_indices), 0,jitter_dist*1.5)
  jitter_angle<-seq(from=angle_from,to=angle_to,length.out=length(label_indices))
  k<-1
  for (i in sorted_label_indices) {
    # Add text label slightly offset from the point
    text_pos<-radial_shift(IMPs1[i],IMPs2[i],jitter_angle[k],dist[k])
    seg_end<-radial_shift(IMPs1[i],IMPs2[i],jitter_angle[k],dist[k]-jitter_dist/2)
    label_base<-strsplit(names(IMPs1)[i],"_")[[1]][length(strsplit(names(IMPs1)[i],"_")[[1]])]
    new_labels <- ifelse(label_base %in% names(replace_map), replace_map[label_base], label_base)
    
    
    text(
      x = text_pos[1],
      y = text_pos[2],
      labels = new_labels,
      cex =1 
    )
    
    # Draw line from point to text label
    segments(
      x0 = IMPs1[i], y0 = IMPs2[i],  # Point coordinates
      x1 = seg_end[1],
      y1 = seg_end[2],
      col = "gray80", lty = 1  
    )
    k=k+1
  }
  
  text(
    x = max(c(IMPs1,IMPs2))*0.13,
    y = max(c(IMPs1,IMPs2))*1,
    labels = paste("rho=",round(cor.test(IMPs1,IMPs2,method="spearman")$estimate,2)),
    cex = 1  # Text size
  )
}

dev.off()
thresholds<-c(80,55,34.3,100)
jitters<-c(30,10,15,40)

# eps 410 x 423
i<-3;RF_correlation_plot_with_label(model_list_rep1[i],model_list_rep2[i],thresholds[i],jitters[i],-80,0)

model_list_rep1_ATX12R37<-list.files(modeldir_rep1,full.names = T)[c(1,2,6,7)]
model_list_Oya2022<-list.files(modeldir_Oya2022,full.names = T)[c(1,3,5,6)]

width_in_inches <- 104.4 / 25.4
height_in_inches <- 110.4 / 25.4
thresholds<-c(75,80,100,100)
jitters<-c(9,12,30,30)


setEPS()
postscript("ATX1.eps", width = width_in_inches, height = height_in_inches, horizontal = FALSE)
i<-4;RF_correlation_plot_with_label(model_list_rep1_ATX12R37[i],model_list_Oya2022[i],thresholds[i],jitters[i],-90,0)
dev.off()
