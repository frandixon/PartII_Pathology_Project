#Plots a heatmap of sample activities

library(ggplot2)

#Import activities
ac <- read.table("activities_proportion_A.txt", header = TRUE, sep = "\t")
ac$SampleFactor <- factor(ac$Samples, levels = ac$Samples)

#Convert numbers to proportions
ac$Proportion_A <- ac$SBS96A/(ac$SBS96A + ac$SBS96B)
ac$Proportion_B <- ac$SBS96B/(ac$SBS96A + ac$SBS96B)

#Convert to plottable format
acP <- data.frame("Sample" = rep(ac$SampleFactor, 2), Signature = c(rep("SBS96A", length(ac[,1])), rep("SBS96B", length(ac[,1]))), Proportion = c(ac$Proportion_A, ac$Proportion_B))

pdf("activity_heatmap_A.pdf")
acPlot <- ggplot(acP, aes(x = Signature, y = Sample, fill = Proportion)) + theme_classic() + geom_tile() + scale_fill_viridis_c() + scale_x_discrete(expand = c(0, 0))
print(acPlot)
dev.off()
