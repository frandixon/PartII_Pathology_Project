#Plots a tanglegram of 2 trees

library(dplyr)
library(ade4)
library(ape)
library(ggtree)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-t1", help = "Tree 1 in newick format")
parser$add_argument("-s", help = "Cosine similarity csv file")
parser$add_argument("-o", help = "Name of output PDF")
args <- parser$parse_args()

#Import the trees and convert to ggtree objects
t1 <- read.tree("genetic_tree_newick.nwk")
#Convert tip labels to uppercase
t1$tip.label <- toupper(t1$tip.label)
p1 <- ggtree(t1) + geom_tiplab(size = 3)
p1 <- flip(p1, 34, 40)
p1 <- flip(p1, 57, 52)
p1 <- flip(p1, 4, 3)
p1 <- flip(p1, 7, 6)
p1 <- flip(p1, 24, 23)
p1 <- flip(p1, 25, 59)
p1 <- flip(p1, 10, 43)
p1 <- flip(p1, 41, 42)
p1 <- flip(p1, 17, 18)

#Import the similarities as a matrix
cosine <- as.matrix(read.csv("cosine_similarity.csv", header = TRUE, row.names = 1))
#Convert labels to uppercase
names(cosine) <- toupper(names(cosine))
row.names(cosine) <- toupper(row.names(cosine))
#Convert the similarities to distances
distances <- 1 - cosine
#Cluster signatures
fit_subs <- hclust(as.dist(distances), method = "average")
#Convert clustering to a tree
fit_subs_tree <- hclust2phylog(fit_subs)
p2 <- ggtree(fit_subs_tree)# + geom_tiplab(size = 4)
p2 <- flip(p2, 48, 45)

#Extract data from ggtree
d1 <- p1$data
d1$x <- d1$x + 0.1
d2 <- p2$data
#Reverse the second tree and move to the right of tree 1
d2$x <- max(d2$x) - d2$x + max(d1$x) + 0.2

#Plot both trees
pp <- p1 + geom_tree(data = d2)

#Add lines between trees
dd <- bind_rows(d1, d2) %>% filter(!is.na(label))
#dd$x <- dd$x + 0.05
#dd$branch <- dd$branch - 0.1
#print(dd[which(dd$label == "CC133_CAT"),])
ppLine <- pp + geom_line(data = dd, aes(x, y, group = label), colour = "grey") #+ geom_text(data = d1, aes(x = x, y = y, label = label), hjust = -0.5) + geom_text(data = d2, aes(x = x, y = y, label = label), hjust = 1)

pdf("tanglegram19.pdf")
print(ppLine)
dev.off()
