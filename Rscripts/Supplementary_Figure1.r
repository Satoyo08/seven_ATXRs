source('custom_functions.r')
#--------  overlap -----------
rep1_bound<-list.files("/Users/Satoyo/Desktop/sevenATXRpaper/seven_ATXRs_localization_git/data/refs/ATX(R)s_bound_gens_list/rep1",full.names = TRUE)
namelist<-paste(c("ATX1","ATX2","ATX3","ATX4","ATX5","ATXR3","ATXR7"),"rep1",sep="_")
for(i in 1:7){
  TMP<-read.table(rep1_bound[i])
  assign(namelist[i], TMP)
}
rep2_bound<-list.files("/Users/Satoyo/Desktop/sevenATXRpaper/seven_ATXRs_localization_git/data/refs/ATX(R)s_bound_gens_list/rep2",full.names = TRUE)
namelist<-paste(c("ATX3","ATX4","ATX5","ATXR3"),"rep2",sep="_")
for(i in 1:4){
  TMP<-read.table(rep2_bound[i])
  assign(namelist[i], TMP)
}
Oya2022_bound<-list.files("/Users/Satoyo/Desktop/temp/論文作業/論文figs/rename_later/Arabidopsis_H3K4me1/data/list_of_genes",
                          full.names = TRUE)[c(1,2,4,5,9,10,8)]
namelist<-c(paste(rep(c("ATX1","ATX2","ATXR7"),each=2),paste(rep(c("rep1","rep2"),3),"Oya22",sep="_"),sep="_")
                  ,"ATXR3_Fiorucci")
for(i in 1:7){
  TMP<-read.table(Oya2022_bound[i])
  assign(namelist[i], TMP)
}

# Inputs for the hyper geometric test
total_population <- 27409
subset_population <- 3000  
drawn_items <- 3000
# lets make euler
library(eulerr)

comparisons <- list(
  list(set1 = ATX3_rep1, set2 = ATX3_rep2, file = "ATX3_rep1_vs_rep2.pdf"),
  list(set1 = ATX4_rep1, set2 = ATX4_rep2, file = "ATX4_rep1_vs_rep2.pdf"),
  list(set1 = ATX5_rep1, set2 = ATX5_rep2, file = "ATX5_rep1_vs_rep2.pdf"),
  list(set1 = ATXR3_rep1, set2 = ATXR3_rep2, file = "ATXR3_rep1_vs_rep2.pdf"),
  list(set1 = ATX1_rep1_Oya22, set2 = ATX1_rep1, file = "ATX1_Oya22_vs_rep1.pdf"),
  list(set1 = ATX2_rep1_Oya22, set2 = ATX2_rep1, file = "ATX2_Oya22_vs_rep1.pdf"),
  list(set1 = ATXR7_rep1_Oya22, set2 = ATXR7_rep1, file = "ATXR7_Oya22_vs_rep1.pdf"),
  list(set1 = ATXR3_Fiorucci, set2 = ATXR7_rep1, file = "ATXR3_Fiorucci_vs_rep1.pdf")
)

compare_list <- function(set1, set2, subset_population=3000, total_population=27409, drawn_items=3000) {
  # Calculate overlap
  overlap <- length(intersect(set1$V1, set2$V1))
  print(overlap)
  # Calculate p-value
  p <- phyper(overlap - 1, subset_population, total_population - subset_population, drawn_items, 
              lower.tail = FALSE, log.p = TRUE)
  
  # Create Euler data
  plot_data <- euler(c(A = length(set1$V1)-overlap, 
                       B = length(set2$V1)-overlap, 
                       "A&B" = overlap))
  # Return overlap and p-value for record
  return(list(overlap = overlap, p_value = p,plotdata=plot_data))
}


i<-1;pdf(comparisons[[i]]$file,width=3,height=2);rs<-compare_list(set1 = comparisons[[i]]$set1,set2 = comparisons[[i]]$set2)
plot(rs$plotdata, quantities = TRUE, labels = c("Rep. 1","Rep. 2"), edges = FALSE,fill = c("#FEE08B", "#3288BD"),alpha = 0.8,main = paste("log10(p-val) =", round(rs$p_value, 2)));dev.off()
i<-2;pdf(comparisons[[i]]$file,width=3,height=2);rs<-compare_list(set1 = comparisons[[i]]$set1,set2 = comparisons[[i]]$set2)
plot(rs$plotdata, quantities = TRUE, labels = c("Rep. 1","Rep. 2"), edges = FALSE,fill = c("#FEE08B", "#3288BD"),alpha = 0.8,main = paste("log10(p-val) =", round(rs$p_value, 2)));dev.off()
i<-3;pdf(comparisons[[i]]$file,width=3,height=2);rs<-compare_list(set1 = comparisons[[i]]$set1,set2 = comparisons[[i]]$set2)
plot(rs$plotdata, quantities = TRUE, labels = c("Rep. 1","Rep. 2"), edges = FALSE,fill = c("#FEE08B", "#3288BD"),alpha = 0.8,main = paste("log10(p-val) =", round(rs$p_value, 2)));dev.off()
i<-4;pdf(comparisons[[i]]$file,width=3,height=2);rs<-compare_list(set1 = comparisons[[i]]$set1,set2 = comparisons[[i]]$set2)
plot(rs$plotdata, quantities = TRUE, labels = c("Rep. 1","Rep. 2"), edges = FALSE,fill = c("#FEE08B", "#3288BD"),alpha = 0.8,main = paste("log10(p-val) =", round(rs$p_value, 2)));dev.off()
i<-5;pdf(comparisons[[i]]$file,width=3,height=2);rs<-compare_list(set1 = comparisons[[i]]$set1,set2 = comparisons[[i]]$set2)
plot(rs$plotdata, quantities = TRUE, labels = c("This study","Oya et al., 2022"), edges = FALSE,fill = c("#FEE08B", "#3288BD"),alpha = 0.8,main = paste("log10(p-val) =", round(rs$p_value, 2)));dev.off()
i<-6;pdf(comparisons[[i]]$file,width=3,height=2);rs<-compare_list(set1 = comparisons[[i]]$set1,set2 = comparisons[[i]]$set2)
plot(rs$plotdata, quantities = TRUE, labels = c("This study","Oya et al., 2022"), edges = FALSE,fill = c("#FEE08B", "#3288BD"),alpha = 0.8,main = paste("log10(p-val) =", round(rs$p_value, 2)));dev.off()
i<-7;pdf(comparisons[[i]]$file,width=3,height=2);rs<-compare_list(set1 = comparisons[[i]]$set1,set2 = comparisons[[i]]$set2)
plot(rs$plotdata, quantities = TRUE, labels = c("This study","Oya et al., 2022"), edges = FALSE,fill = c("#FEE08B", "#3288BD"),alpha = 0.8,main = paste("log10(p-val) =", round(rs$p_value, 2)));dev.off()
i<-8;pdf(comparisons[[i]]$file,width=3,height=2);rs<-compare_list(set1 = comparisons[[i]]$set1,set2 = comparisons[[i]]$set2)
plot(rs$plotdata, quantities = TRUE, labels = c("This study","Fiorucci et al., 2022"), edges = FALSE,fill = c("#FEE08B", "#3288BD"),alpha = 0.8,main = paste("log10(p-val) =", round(rs$p_value, 2)));dev.off()

