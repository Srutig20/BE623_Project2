getwd()
setwd("C:/Users/USER/Downloads/project")
#loading the libraries/packages
library("DESeq2")
library("dplyr")
library("ggplot2")
library("pheatmap")

#TASK 1a

#reading the csv files
read.csv("airway_counts.csv")
counts <- read.csv("airway_counts.csv", row.names = 1)
read.csv("airway_samples.csv")
samples <- read.csv("airway_samples.csv")
View(counts)

#reporting dimensions of counts
dim(counts)

#showing structure of samples
str(samples)

#TASK 1b
#just mutating does not mutate the original file
mutate(samples, lib_size=colSums(counts))
#assigning the mutation to the samples data again adds the
#change to the original file
samples <- mutate(samples, lib_size=colSums(counts))
#sample table can be viewed
View(samples)
# to compute per group statistics
samples %>%
  group_by(group)%>%
  summarize(
    mean_library_size = mean(lib_size),
    n = n()
  )


#TASK 1c

ggplot(samples, aes(x=cell, y =lib_size, fill=group)) + geom_col()+
  labs(
    title = "Library size of cells",
    x = "cells",
    y = "Library_size"
  )

ggsave("plots/library_sizes.png", width = 6, height = 4, dpi = 300)

#TASK 2a

rownames(samples) <- samples$run
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                              colData = samples, design = ~ group)

levels(dds$group)

#TASK 2b
#removing low count genes
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
nrow(dds)

#TASK 2c
vsd <- vst(dds, blind = TRUE)
dir.create("plots", showWarnings = FALSE)
p_pca <- plotPCA(vsd, intgroup = "group")
p_pca

ggsave("plots/PCA.png", p_pca, width = 7, height = 5, dpi = 300) 
#%variance of PC1 = 42%, %variance of PC2 = 26%
# Yes the trt and untrt samples separate along x-axis.

#TASK 3a

dds <- DESeq(dds)
# Normalised counts (useful to look at, or to plot single genes)
write.csv(counts(dds, normalized = TRUE), "results/normalized_counts.csv", quote = FALSE)

#TASK 3b
comparison <- paste0("trt", "_vs_", "untrt")
res <- results(dds, contrast = c("group", "trt", "untrt"))
head(res)               

res <- res[order(res$padj), ]
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$status <- "Not significant"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- "Up"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "Down"

table(res_df$status)
#up = 465
#down = 388

degs <- subset(res_df, status %in% c("Up", "Down"))
nrow(degs)

#TASK 3c
write.csv(res_df, paste0("results/", comparison, "_all_genes.csv"), row.names = FALSE)
write.csv(degs, paste0("results/", comparison, "_DEGs.csv"), row.names = FALSE)
deg_table <- read.csv(paste0("results/", comparison, "_DEGs.csv"))
head(deg_table)
head(deg_table, 5)

#TASK 3d
#i
# Answer: We give DESeq2 raw counts and not normalized counts because it does the
# normalisation itself.
#ii
#Answer: 
#iii
# Answer: If a gene has log2FC =3, it means that the fold value is 2^3 = 8.

#TASK 4a
volcano_data <- res_df[!is.na(res_df$padj), ] # remove genes with padj = NA
p_volcano <- ggplot(volcano_data,
                    aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up" = "firebrick",
                                "Down" = "steelblue",
                                "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = comparison,
       x = "log2 fold change",
       y = "-log10 adjusted p-value") +
  theme_bw()
p_volcano 
#Interpretation: The plot shows downregulation and upregulation of treated vs 
#untreated cells. The blue dots above the horizontal line(padj 0.05; significance level)
#are considered to be statistically significant data points and are downregulated.
#The grey points are considered to be not significants points and the red points
#are considered to be upregulated by the treatment.

# TASK 4b
str(counts)
counts[1=="ENSG00000120129"]
