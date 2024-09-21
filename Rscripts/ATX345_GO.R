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

