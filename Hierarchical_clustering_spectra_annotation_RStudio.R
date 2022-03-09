#Calculates a plots hierarchical clustering of a set of cosine similarities between spectra

library("gplots")
library("ade4")
library("ggtree")
library("ggplot2")

#Import the similarities as a matrix
cosine <- as.matrix(read.csv("cosine_similarity.csv", header = TRUE, row.names = 1))
#Convert the similarities to distances
distances <- 1 - cosine

#Import annotations to plot next to the tree
annotations <- read.csv("sample_annotations.csv", header = TRUE, row.names = 1)

#Set colours for the annotations
colours <- read.csv("annotation_colours.csv", header = TRUE)
labels <- list()
for (i in 1:length(colours[,1])) {
	labels[colours[i,1]] <- colours[i,2]}

outFile <- "hierarchical_clustering"

#Cluster signatures
fit_subs <- hclust(as.dist(distances), method = "average")

#Convert clustering to a tree
fit_subs_tree <- hclust2phylog(fit_subs)

pdf(paste(outFile, ".pdf", sep = ""))
p <- ggtree(fit_subs_tree) + geom_tiplab()
tPlot <- gheatmap(p, annotations, offset = 0.05, width = 0.25, font.size = 0) + scale_fill_manual(values = labels) + theme(legend.title = element_blank())
print(tPlot)
dev.off()
