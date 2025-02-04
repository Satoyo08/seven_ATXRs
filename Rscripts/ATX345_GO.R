# Load necessary libraries
library(ggplot2)
library(scales)

# Read 
Rdata <- read.table("../data/ATX345-marked_GO.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE,fill = TRUE,skip=12, quote = "")
head(data)
colnames(Rdata)<-c("BiologicalProcess","inputlist","total","expected","over/under","Fold.Enrichment","Raw.P.value")

Rdata$FoldEnrichment <- as.numeric(Rdata$Fold.Enrichment) 
Rdata$PValue <- as.numeric(Rdata$Raw.P.value) 
nsdata<-Rdata[Rdata$Fold.Enrichment>2.5,]
data<-nsdata[rev(order(nsdata$PValue)),]

# BiologicalProcess sorted by PValue and make it factor
data$BiologicalProcess <- factor(data$BiologicalProcess, levels = unique(data$BiologicalProcess))


# Create the plot
p <- ggplot(data, aes(x = "", y = BiologicalProcess)) +  # Blank x-axis
  geom_point(aes(size = FoldEnrichment, color = PValue)) + 
  scale_size(range = c(3, 15)) +  # Adjust bubble size
  scale_color_gradient(low = "red", high = "blue", trans = "log10", labels = scales::scientific) +  # Color based on p-value
  theme_classic() +
  labs(
    title = "ATX3/4/5/-marked genes Ontology enrichment",
    x = NULL,
    y = "Biological Process",
    color = "P-value (log scale)",
    size = "Fold Enrichment"
  ) +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
    axis.title.x = element_blank()  # Hide x-axis title
  )

# Print the plot
print(p)




### GO for Chen 2017 dataset

library(readxl)
class1 <- read_excel("/Users/Satoyo/Downloads/plphys_v174_3_1795_s1 (2)/PP2016-01944R2_Supplemental_Tables_.xlsx", sheet = "Table S6")
class2 <- read_excel("/Users/Satoyo/Downloads/plphys_v174_3_1795_s1 (2)/PP2016-01944R2_Supplemental_Tables_.xlsx", sheet = "Table S7")

# rm header
class1_dat <- class1[-c(1, 2),c(1:3) ]; colnames(class1_dat)<-c('chr','st','en')
class2_dat <- class2[-c(1, 2), ];colnames(class2_dat)<-c('chr','st','en')

# merge
data<-rbind(class1_dat,class2_dat)

write.table(data, file = "../data/refshypoH3K4me2_sites_from_Chen2017.bed", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
system('bedtools intersect -wa -a /Users/Satoyo/araport11_beds/araport11_all_sorted_nuclear.bed -b hypoH3K4me2_sites_from_Chen2017.bed | cut -f 4 > Chen_hypoH3K4me2_overlapgenes_AtID.txt')

# do GO in panther

Rdata <- read.table("../data/refs/hypoH3K4me2_gene_from_Chen2017_GO.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE,fill = TRUE,skip=12, quote = "")
head(Rdata)
colnames(Rdata)<-c("BiologicalProcess","inputlist","total","expected","over/under","Fold.Enrichment","Raw.P.value")

Rdata$FoldEnrichment <- as.numeric(Rdata$Fold.Enrichment) 
Rdata$PValue <- as.numeric(Rdata$Raw.P.value) 
nsdata<-Rdata[Rdata$Fold.Enrichment>2.5,]
data<-nsdata[rev(order(nsdata$PValue)),]

# BiologicalProcess sorted by PValue and make it factor
data$BiologicalProcess <- factor(data$BiologicalProcess, levels = unique(data$BiologicalProcess))


# Create the plot
p <- ggplot(data, aes(x = "", y = BiologicalProcess)) +  # Blank x-axis
  geom_point(aes(size = FoldEnrichment, color = PValue)) + 
  scale_size(range = c(3, 15)) +  # Adjust bubble size
  scale_color_gradient(low = "red", high = "blue", trans = "log10", labels = scales::scientific) +  # Color based on p-value
  theme_classic() +
  labs(
    title = "ATX3/4/5/-marked genes Ontology enrichment",
    x = NULL,
    y = "Biological Process",
    color = "P-value (log scale)",
    size = "Fold Enrichment"
  ) +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
    axis.title.x = element_blank()  # Hide x-axis title
  )

# Print the plot
print(p)

# also enriched in hormone metabolics
