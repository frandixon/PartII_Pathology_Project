#Compares pairs of spectra to identify significantly different mutations

library(tidyverse)
library(data.table)
library(ape)
library(ggthemes)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-s", nargs = "+", help = "Spectra to be compared, all pairs will be compared")
parser$add_argument("-o", help = "Prefix of output plots")
args <- parser$parse_args()

fileList <- c("file1.csv", "file2.csv")

cnames <- c("Sample", "Branch", "Branch_length", "Total_mutations", "C>A", "C>G", "C>T", 
            "T>A", "T>C", "T>G", "C>A_proportion", "C>G_proportion", "C>T_proportion", "T>A_proportion", 
            "T>C_proportion", "T>G_proportion")

#Import the branches into a single data frame
data <- map_dfr(fileList, ~{
  f <- fread(.x, skip = 1, drop = 16) %>% as_tibble() %>% add_column(sample = gsub(".*/", 
                                                                                 "", gsub("_branches.*", "", .x)), .before = 1)
  colnames(f) <- cnames
  return(f)})

temp_data <- data

#6 mutation types
muts <- c("`C>A`", "`C>G`", "`C>T`", "`T>A`", "`T>C`", "`T>G`")

#Compare the regression between each mutation and branch length between each pair of samples
results <- map_dfr(muts, ~{
  #print(.x)
  temp_data$Sample <- as.factor(temp_data$Sample)
  form <- as.formula(paste(.x, "Sample + offset(log(Total_mutations))", sep = " ~ "))
  # form <- as.formula(paste(.x, 'GPSC + Branch_length +
  # offset(log(Total_mutations))', sep = ' ~ '))
  
  #Fit linear model between mutation and branch length
  m <- glm(form, data = temp_data, family = "poisson")
  # m <- MASS::glm.nb(form, data = temp_data)
  
  s <- summary(multcomp::glht(m, multcomp::mcp(Sample = "Tukey")))
  
  tibble(mut = .x, pair = names(s$test$coefficients), coefficients = s$test$coefficients, 
         `Std. Error` = s$test$sigma, `z value` = s$test$tstat, p.value = s$test$pvalues)
  
}) %>% arrange(p.value)

results$adj.p.value <- p.adjust(results$p.value, method = "BH")
results$mut <- gsub("`", "", results$mut)
tophits <- results %>% filter(adj.p.value < 0.05)
tophits$p.value <- format(tophits$p.value, digits = 3)
tophits$adj.p.value <- format(tophits$adj.p.value, digits = 3)
knitr::kable(tophits)

#Plot proportion of each mutation type in each sample
pdf(paste("prefix", "_mutation_proportions.pdf", sep = ""))
pdf <- data[, c(1:4, 11:ncol(data))] %>% tidyr::pivot_longer(cols = colnames(data)[11:ncol(data)])
pdf$mutation <- gsub("_proportion", "", pdf$name)

pdf <- pdf %>% group_by(Sample, mutation) %>% summarise(value = mean(value))

pPlot <- ggplot(pdf, aes(x = mutation, y = value, fill = mutation)) + geom_col() + facet_wrap(~Sample) + 
  theme_clean(base_size = 20) + theme(plot.background = element_blank(), legend.background = element_blank(), 
                                      axis.title.x = element_blank()) + ylab("Proportion of mutations")
print(pPlot)
dev.off()

#Plot differences in means of each mutation type
pdf(paste("prefix", "_difference_in_means.pdf", sep = ""))
pdf2 <- results
sp <- apply(str_split_fixed(pdf2$pair, " - ", 2), 1, sort)
levs <- sort(unique(c(sp)))
pdf2$SampleA <- factor(sp[1, ], levels = levs)
pdf2$SampleB <- factor(sp[2, ], levels = levs)
pdf2$diff <- pdf$value[match(paste(pdf2$mut, pdf2$SampleA), paste(pdf$mutation, pdf$Sample))] - 
  pdf$value[match(paste(pdf2$mut, pdf2$SampleB), paste(pdf$mutation, pdf$Sample))]

mPlot <- ggplot(pdf2, aes(x = mut, y = diff, fill = mut, alpha = adj.p.value < 0.05)) + geom_col() + 
  facet_grid(SampleA ~ SampleB, drop = TRUE) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 90, 
                                                                                                       vjust = 0.2, hjust = 1)) + scale_alpha_manual(values = c(0.2, 1)) + ylab("difference in mean") + 
  xlab("mutation")
print(mPlot)
dev.off()

muts <- c("`C>A`", "`C>G`", "`C>T`", "`T>A`", "`T>C`", "`T>G`")

results <- map_dfr(unique(data$Sample), function(samples) {
  temp_data <- data
  temp_data$Sample[temp_data$Sample != samples] <- "all"
  res <- map_dfr(muts, ~{
    form <- as.formula(paste(.x, "Sample + Branch_length + offset(log(Total_mutations))", 
                             sep = " ~ "))
    broom::tidy(MASS::glm.nb(form, data = temp_data)) %>% filter(grepl("Sample", 
                                                                       term)) %>% add_column(mut = .x, .before = 1)
    # broom::tidy(glm(form, data = temp_data, family = 'quasipoisson')) %>%
    # filter(grepl('GPSC', term)) %>% add_column(mut=.x, .before = 1)
  }) %>% add_column(Sample = samples, .before = 1)
  return(res)
}) %>% arrange(p.value)

results$mut <- gsub("`", "", results$mut)
results$adj.p.value <- p.adjust(results$p.value, method = "BH")

tophits <- results %>% filter(adj.p.value < 0.05)
tophits$p.value <- format(tophits$p.value, digits = 3)
tophits$adj.p.value <- format(tophits$adj.p.value, digits = 3)
knitr::kable(tophits)
