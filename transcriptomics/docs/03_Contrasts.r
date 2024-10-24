library(eulerr)
library(DESeq2)
library(ggplot2)

# picking up from import of counts matrix and meta data and creation of dds object from previous script:
# "Transcriptomics_Days1-4.R"

options(bitmapType="cairo")

# Import the counts matrix
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt", header=TRUE, row.names=1)
countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)


dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ DevTemp + FinalTemp)
dds <- dds[rowSums(counts(dds) >= 10) >= 15,]

# Let's set up references for each contrast between treatments 

dds$treatment <- relevel(dds$DevTemp, ref = "D18")
dds$generation <- relevel(dds$FinalTemp, ref = "BASE")


# set reference levels
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))

design(dds) <- ~group
dds <- DESeq(dds)

dim(dds)
resultsNames(dds)
# 
# [1] "Intercept"               "group_D18A33_vs_D18A28"  "group_D18BASE_vs_D18A28"
# [4] "group_D22A28_vs_D18A28"  "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

#res_D18_BASE_A28 <- results(dds, contrast=c("group", "D18BASE", "D18A28"), alpha=0.05)

# 1. compare baseline gene expression between developmental treatment groups
res_D18_BASE_D22_BASE <- results(dds, contrast=c("group", "D18BASE", "D22BASE"), alpha=0.05)
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),]
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),]
head(res_D18_BASE_D22_BASE)
summary(res_D18_BASE_D22_BASE)

# make a list of differentially expressed genes (DEGs) and plot
degs_D18_BASE_D22_BASE <- row.names(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj<0.05,])
plotMA(res_D18_BASE_D22_BASE, ylim=c(-4,4))

# 2. Now compare gene expression between developmental treatment groups @ A28
res_D18_A28_D22_A28 <- results(dds, contrast=c("group", "D18A28", "D22A28"), alpha=0.05)
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[!is.na(res_D18_A28_D22_A28$padj),]
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[order(res_D18_A28_D22_A28$padj),]
head(res_D18_A28_D22_A28)
summary(res_D18_A28_D22_A28)

# make a list of differentially expressed genes (DEGs) and plot
degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj<0.05,])
plotMA(res_D18_A28_D22_A28, ylim=c(-4,4))

# 3. Now compare gene expression between developmental treatment groups @ A33
res_D18_A33_D22_A33 <- results(dds, contrast=c("group", "D18A33", "D22A33"), alpha=0.05)
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[!is.na(res_D18_A33_D22_A33$padj),]
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[order(res_D18_A33_D22_A33$padj),]
head(res_D18_A33_D22_A33)
summary(res_D18_A33_D22_A33)

# make a list of differentially expressed genes (DEGs) and plot
degs_D18_A33_D22_A33 <- row.names(res_D18_A33_D22_A33[res_D18_A33_D22_A33$padj<0.05,])
plotMA(res_D18_A33_D22_A33, ylim=c(-4,4))

length(degs_D18_BASE_D22_BASE)
length(degs_D18_A28_D22_A28)
length(degs_D18_A33_D22_A33)

length(intersect(degs_D18_BASE_D22_BASE, degs_D18_A28_D22_A28)) #107

length(intersect(degs_D18_BASE_D22_BASE, degs_D18_A33_D22_A33)) #44

length(intersect(degs_D18_A33_D22_A33, degs_D18_A28_D22_A28)) #29

int <- intersect(degs_D18_BASE_D22_BASE, degs_D18_A28_D22_A28)
length(intersect(degs_D18_A33_D22_A33, int)) #23

1935-107-44+23 # 1807 genes diff expressed uniquely at baseline btwn 18 vs. 22
296-107-29+23 # 183 uniquely expressed when exposed to A28
778-44-29+23 # 28 uniquely expressed when exposed to A33

107-23 # 84 unique to baseline and A28
44-23 # 21 genes unique to baseline and A33
29-23 # 6 unique to A28 and A33

# Now have all unique pieces needed to make Euler plot:

myEuler <- euler(c("BASE"=1807, 
                   "A28"=28, 
                   "BASE&A28"=84, 
                   "BASE&A33"=21,
                   "A28&A33"=6,
                   "BASE&A28&A33"=23))

plot(myEuler, lty=1:3, quantities=T)

########################################

# contrast D18_A28vsBASE

res_D18_BASEvsA28 <- as.data.frame(results(dds, 
                                           contrast=c("group","D18BASE","D18A28"), alpha=0.05))
# contrast D22_A28vsBASE

res_D22_BASEvsA28 <- as.data.frame(results(dds, 
                                           contrast=c("group","D22BASE","D22A28"), alpha=0.05))

# merge dataframes

res_df28 <- merge(res_D18_BASEvsA28, res_D22_BASEvsA28, by="row.names", suffixes=c(".18",".22"))

rownames(res_df28) <- res_df28$Row.names

res_df28 <- res_df28[,-1]

library(tidyverse)

# Define colors with mutate function

res_df28 <- res_df28 %>% 
  mutate(fill=case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "turquoise2",
    padj.18 < 0.05 & stat.18 > 0 ~ "magenta1",
    padj.22 < 0.05 & stat.18 < 0 ~ "blue2",
    padj.22 < 0.05 & stat.18 > 0 ~ "red"
  ))

color_counts <- res_df28 %>%
  group_by(fill) %>%
  summarise(count = n())

label_positions <- data.frame(fill=c("blue2","magenta1","red","turquoise2"),
                              x_pos=c(1,5,0,-7.5),
                              y_pos=c(-5,0,9,3))


label_data <- merge(color_counts,label_positions,by="fill")

# plot!

ggplot(res_df28, aes(x=log2FoldChange.18, y=log2FoldChange.22, color=fill)) +
  geom_point(alpha=1) +
  scale_color_identity() +
  geom_text(data=label_data, aes(x=x_pos, y=y_pos, label=count, color=fill), size=5) +
  geom_abline(intercept=0, slope=1, linetype="dashed",color="black") +
  geom_abline(intercept=0, slope=-1, linetype="dashed",color="gray") +
  xlim(-10,10) + ylim(-10,10) +
  labs(x="Log2FoldChange 28 vs. BASE at 18",
       y="log2FoldChange 28 vs. BASE at 22", title="How does response to 28C vary by DevTemp?") +
  theme_minimal()

###################################################
# NOT FINISHED MODIFYING CODE BELOW YET...NEED TO CHG TO A33

# Now do same for comparison at A33

# contrast D18_A33vsBASE

res_D18_BASEvsA33 <- as.data.frame(results(dds, 
                                           contrast=c("group","D18BASE","D18A33"), alpha=0.05))
# contrast D22_A33vsBASE

res_D22_BASEvsA33 <- as.data.frame(results(dds, 
                                           contrast=c("group","D22BASE","D22A33"), alpha=0.05))

# merge dataframes

res_df33 <- merge(res_D18_BASEvsA33, res_D22_BASEvsA33, by="row.names", suffixes=c(".18",".22"))

rownames(res_df33) <- res_df33$Row.names

res_df33 <- res_df33[,-1]

# Define colors with mutate function

res_df33 <- res_df33 %>% 
  mutate(fill=case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "turquoise2",
    padj.18 < 0.05 & stat.18 > 0 ~ "magenta1",
    padj.22 < 0.05 & stat.18 < 0 ~ "blue2",
    padj.22 < 0.05 & stat.18 > 0 ~ "red"
  ))

# plot!

ggplot(res_df33, aes(x=log2FoldChange.18, y=log2FoldChange.22, color=fill)) +
  geom_point(alpha=1) +
  scale_color_identity() +
  geom_text(data=label_data, aes(x=x_pos, y=y_pos, label=count, color=fill), size=5) +
  geom_abline(intercept=0, slope=1, linetype="dashed",color="black") +
  geom_abline(intercept=0, slope=-1, linetype="dashed",color="gray") +
  xlim(-10,10) + ylim(-10,10) +
  labs(x="Log2FoldChange 33 vs. BASE at 18",
       y="log2FoldChange 33 vs. BASE at 22", title="How does response to 33C vary by DevTemp?") +
  theme_minimal()






