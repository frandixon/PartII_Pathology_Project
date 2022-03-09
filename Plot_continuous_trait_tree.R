#Plots a continuous trait next to a tree

library("ape")
library("ggtree")
library("ggplot2")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("-t", help = "Tree next to which trait will be plotted")
parser$add_argument("-a", help = "Annotations to be plotted next to the tree, first column contains sample names to match those in the tree")
parser$add_argument("-o", help = "Name of output pdf")
args <- parser$parse_args()

#Import the tree and convert to ggtree object
t <- read.tree("genetic_tree_newick.nwk")
p <- ggtree(t)

#Import trait to be plotted
trait <- read.csv("labels_file.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

pdf("output.pdf")
tPlot <- gheatmap(p, trait, offset = 0.08, width = 0.75) + geom_tiplab() + scale_fill_viridis_c()
print(tPlot)
dev.off()
