#Clusters signatures based on their cosine distance
#Takes a table of cosine similarities
#Calculates distances as 1 - cosine similarity

library(argparse)

plotPartitions <- function(hc_fit, partitions, cutoff, outFile) {
  #Number of signatures in each partition
  table_n_part <- table(partitions)
  #Reorder table by partition size
  #new_partitions_order <- order(table_n_part, decreasing = TRUE)
  
  #partitions_copy <- partitions
  #for (i in 1:length(partitions_copy)) {
  #  partitions[partitions_copy == new_partitions_order[i]] <- i}
  
  partition_list <- unique(partitions)
  boxes_order <- partitions[hc_fit$order]
  
  jpeg(outFile, width = 2000, height = 1000, res = 160)
  plot(fit_subs, hang = -1,xlab = "",ylab = "1 - cosine similarity, average linkage",sub = "")
  abline(a = cutoff, b = 0)
  start_draw <- 0.5
  box_bottom <- -0.62
  for (i in unique(boxes_order)) {
    end_draw <- start_draw + sum(boxes_order == i)
    par(xpd = TRUE)
    rect(start_draw, box_bottom, end_draw, 0, border = "red")
    text(i, x = start_draw + (end_draw - start_draw)/2, y = box_bottom - 0.05, cex = 0.8)
    par(xpd = FALSE)
    start_draw <- end_draw
  }
  dev.off()

}
