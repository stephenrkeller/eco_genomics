### 01_DESeq2.R
### Script for importing our counts data, exploring the data, and testing for differential expression.

setwd("~/projects/eco_genomics/transcriptomics") # for in class on server
#setwd("/Users/mpespeni/Dropbox/2_Teaching/Intro_Ecological_Genomics/from_Alison") # when on my machine

library(DESeq2)
library(ggplot2)

options(bitmapType="cairo")

# Import the counts matrix
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt", header=TRUE, row.names=1)
#countsTable <- read.table("tonsa_counts.txt", header=TRUE, row.names=1)

head(countsTable)
dim(countsTable) #119438 21

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
#added letters to the number values ahead of the DevTemp and AcuteTemp values in experimental details. This is because we got the error that the matrix is unbalanced, we think because the dev temp and acute temp have the same vales sometimes. 

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

# Let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound)) #18454529

barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound),cex.names=0.5, las=3,ylim=c(0,20000000))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd=2)

# the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # 3244.739 # [1] 11930.81 - tonsa, 6076.078 - hudsonica genes, 2269 - hudsonica isoform
median(rowSums(countsTableRound)) # 64 # [1] 2226 - tonsa, 582 - hudsonica, 109

apply(countsTableRound,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound,1,mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound,1,mean),xlim=c(0,1000), ylim=c(0,120000),breaks=10000)

#### Create a DESeq object and define the experimental design here with the tilda
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ DevTemp + FinalTemp)


dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 10 in more than 75% of samples, so ~15)
## suggested by WGCNA on RNAseq FAQ


dds <- dds[rowSums(counts(dds) >= 10) >= 15,] #15 is 75% of 21
nrow(dds) #35527 number of transcripts with more than 10 reads in more than or equal to 15 samples

# Run the DESeq model to test for global differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)



# PCA to visualize global gene expression patterns
# first transform the data for plotting using variance stabilization
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("DevTemp","FinalTemp"), returnData=TRUE)
percentVar <- round(100*attr(pcaData,"percentVar")) 
final_temp_colors <- c("BASE" = "darkolivegreen4", "A28" = "darkorange2", "A33" = "darkred")
shapes_choose <- c("D18" = 16, "D22" = 18)

p <- ggplot(pcaData, aes(PC1, PC2, color = FinalTemp, shape = DevTemp)) +
  geom_point(size = 5) +  # Increased dot size
  scale_shape_manual(values = shapes_choose) +
  scale_color_manual(values = final_temp_colors) +
  labs(x = paste0('PC1: ', percentVar[1], ' %'),
       y = paste0('PC2: ', percentVar[2], ' %')) +
  theme_bw(base_size = 16) +  # Increased base font size
  coord_fixed() +
  theme(
    panel.border = element_rect(colour = "black", size = 1.5),  # Bold outer box
    axis.text = element_text(size = 14),  # Font size for axis text
    axis.title = element_text(size = 16),  # Font size for axis titles
    legend.title = element_text(size = 14),  # Font size for legend titles
    legend.text = element_text(size = 12)  # Font size for legend text
  ) +
  guides(
    shape = guide_legend(title = "DevTemp", reverse = FALSE, order = 1),
    color = guide_legend(title = "FinalTemp", reverse = FALSE, order = 2)
  )

p


####################################################
#
# Day 3 of transcriptomics - Test for differential expression & visualize
#
####################################################
library(pheatmap)
# Recall the results you've generated
resultsNames(dds)
# [1] "Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28" 
# [4] "FinalTemp_BASE_vs_A28"

# pull out the results for D22 vs D18 and set alpha = 0.05
res_D22vsD18 <- results(dds, name="DevTemp_D22_vs_D18", alpha=0.05) 

# order by significance
res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]
head(res_D22vsD18) 
# log2 fold change (MLE): DevTemp D22 vs D18 
# Wald test p-value: DevTemp D22 vs D18 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN140854_c0_g5_i2  1210.982      -1.195547 0.0783055 -15.26772 1.25475e-52 4.28270e-48
# TRINITY_DN142149_c1_g3_i1   416.434      -1.648860 0.1200891 -13.73031 6.68496e-43 1.14086e-38
# TRINITY_DN140616_c0_g2_i1   835.221      -0.956038 0.0833345 -11.47229 1.81787e-30 2.06825e-26
# TRINITY_DN140854_c0_g5_i1  1382.932       0.818620 0.0769785  10.63439 2.06178e-26 1.75932e-22
# TRINITY_DN147593_c1_g4_i1   103.376       1.410794 0.1583056   8.91184 5.01926e-19 3.42635e-15
# TRINITY_DN135172_c0_g1_i3   129.693       1.732298 0.2090688   8.28578 1.17336e-16 6.67486e-13

# NOTE: The orientation of the levels in the first line indicates the direction 
# of the comparison; D22 vs D18 means that positive log2FoldChange and stat values 
# have higher expression in D22 than D18, negative log2FoldChange and stat values
# have higher expression in D18 than D22.

summary(res_D22vsD18)
# out of 35527 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2249, 6.3%
# LFC < 0 (down)     : 2235, 6.3%
# outliers [1]       : 17, 0.048%
# low counts [2]     : 1378, 3.9%

# Note the set adjusted p-value is correct.

### Plot Individual genes ### 

# Counts of specific top gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="TRINITY_DN140854_c0_g5_i2", intgroup = (c("DevTemp","FinalTemp")), returnData=TRUE)
d

p <-ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape=FinalTemp)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p



### We can make an MA plot to see the relationship between LFC and counts

plotMA(res_D22vsD18, ylim=c(-4,4))


### We can make a volcano plot to see the relationship between LFC and p-value

# Convert results to a data frame for easy plotting
res_df <- as.data.frame(res_D22vsD18)

# Add a column to label significantly differentially expressed genes
# Adjust the thresholds for significance (e.g., log2 fold change > 1 and adjusted p-value < 0.05)
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")

# Make the volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
  geom_point(alpha = 0.8) +  # Scatter plot with transparency
  scale_color_manual(values = c("grey", "purple")) +  # Customize colors
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot") +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Horizontal line for p-value threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue")  # Vertical lines for fold change threshold



### We can also make a heatmap of the DEGs

# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)

vsd <- vst(dds, blind=FALSE)

topgenes <- head(rownames(res_D22vsD18),20) # pull out the top 20 genes
mat <- assay(vsd)[topgenes,] # make a matrix of the data for these genes across all samples
mat <- mat - rowMeans(mat) # normalize across genes by subtracting the mean 
df <- as.data.frame(colData(dds)[,c("DevTemp","FinalTemp")])
pheatmap(mat, annotation_col=df, show_rownames=FALSE, cluster_cols=T,cluster_rows = T)

# NOTE: you should see clustering by devtemp because these are the top devtemp genes
# What if you did all of the above for the "FinalTemp_BASE_vs_A28" contrast?


# Let's set up references for each contrast between treatments 
. 
dds$treatment <- relevel(dds$DevTemp, ref = "D18")
dds$generation <- relevel(dds$FinalTemp, ref = "BASE")


dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group
dds <- DESeq(dds)
dim(dds)
resultsNames(dds)
# [1] "Intercept"               "group_D18A33_vs_D18A28"  "group_D18BASE_vs_D18A28"
# [4] "group_D22A28_vs_D18A28"  "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

#contrast D18_A28
res_D18_BASE_A28 <- results(dds, contrast=c("group","D18BASE","D18A28"), alpha = 0.05)

res_D18_BASE_A28 <- res_D18_BASE_A28[order(res_D18_BASE_A28$padj),]
head(res_D18_BASE_A28)
summary(res_D18_BASE_A28)
# LFC > 0 (up)       : 11, 0.031%
# LFC < 0 (down)     : 30, 0.084%
# outliers [1]       : 3, 0.0084%
# low counts [2]     : 0, 0%

# Let's see what a top gene looks like - if it's what we expect

d <-plotCounts(dds, gene="TRINITY_DN70277_c0_g1_i1", intgroup = (c("DevTemp","FinalTemp")), returnData=TRUE)
d

p <-ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape=FinalTemp)) +
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p

# Plot again zoomed in... Now we can see the difference between BASE and A28
p <-ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape=FinalTemp)) + ylim(0,1100) +
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p


# Can continue to make series of euler plots and a panel of volcano plots of all the contrasts
