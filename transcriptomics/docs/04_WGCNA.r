### Script for analyzing and visualizing gene coexpression networks

library(ggplot2)
library(DESeq2)
library(WGCNA); options(stringsAsFactors=F)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

options(bitmapType="cairo")

setwd("~/projects/eco_genomics/transcriptomics/")

# Step 1: Import the data
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt", header=TRUE, row.names=1)
countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)

# Import trait data
traitData <- read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", header=T, row.names=1)

#Filter the matrix to just BASE data (b/c that's what we have traits for)
filtered_count_matrix_BASEonly <- countsTable[,conds$FinalTemp == "BASE"]
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp == "BASE",]
rounded_filtered_count_matrix <- round(filtered_count_matrix_BASEonly)

# Step 2: detecting outliers
# detect outlier genes
gsg <- goodSamplesGenes(t(rounded_filtered_count_matrix))
summary(gsg)

table(gsg$goodGenes)
# FALSE  TRUE 
# 37235 82203 

table(gsg$goodSamples) # All good!

data_WGCNA <- rounded_filtered_count_matrix[gsg$goodGenes==T,]
dim(data_WGCNA)

# use clustering with a tree dendrogram to identify outlier samples
htree <- hclust(dist(t(data_WGCNA)), method="average")
plot(htree)

# PCA outlier detection
pca <- prcomp(t(data_WGCNA))
pca_data <- data.frame(pca$x)

pca.var <- pca$sdev^2
pca.var.pct <- round(pca.var/sum(pca.var)*100,2)

ggplot(pca_data, aes(x=PC1, y=PC2)) +
  geom_point() +
  geom_text(label=rownames(pca_data)) +
  labs(x=paste0("PC1: ",pca.var.pct[1]," %"),
       y=paste0("PC2: ",pca.var.pct[2]," %"))

# based on the initial tests, the samples seemed fine so we'll leave them in (?)...

# Step 3:
# Now, want deseq to cluster without any info about groups...so design an intercept only model

dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA,
                                    colData = filtered_sample_metadata_BASEonly,
                                    design= ~1) # there are no specified groups

dds_WGCNA_75 <- dds_WGCNA[rowSums(counts(dds_WGCNA)>=15) >=6, ] # sum of counts for each gene has to be at least 15 for at least 6 samples
nrow(dds_WGCNA_75)

dds_norm <- vst(dds_WGCNA_75) # perform variance stablization

# get and save normalized counts to use below
norm.counts <- assay(dds_norm) %>%
  t()

# Step 4: Network construction

# choose a set of soft-thresholding powers

power <- c(c(1:10), seq(12,50,2))

# Call the network topology analysis fcuntion (takes a few mins to run)

sft <- pickSoftThreshold(norm.counts, 
                         powerVector = power,
                         networkType = "signed",
                         verbose=5) 

# signed means we're only paying attention to genes that are positively correlated 
# with each other -- either both positive or both negative logFC

sft.data <- sft$fitIndices

# plot to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label=Power)) +
  geom_point() +
  geom_text(nudge_y=0.1) +
  geom_hline(yintercept=0.8, color="red") +
  labs(x="Power",y="Scale free topology model fit, signed R^2") +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label=Power)) +
  geom_point() +
  geom_text(nudge_y=0.1) +
  geom_hline(yintercept=0.8, color="red") +
  labs(x="Power",y="Mean Connectivity") +
  theme_classic()

grid.arrange(a1,a2,nrow=2)

# WGCNA wants us to pick a relationship across the many genes and 6 samples
# as you increase power, you increase the model fit of the toplpgy
# Want to choose a threshold of power that has high fit but need to be careful about
# the decrease in connectivity...its a tradeoff. Here, a value around 24 looks good...



