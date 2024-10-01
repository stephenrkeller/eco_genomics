# This script will test for selection outliers using the PCAdapt method
# this is an extension of the Fst outlier method, developed by Duforest-Fribourg et al. (2015)
# with updates by Luu et al. (2017) and Prive et al. (2020)
# Github page here: https://github.com/bcm-uga/pcadapt

library(tidyverse)
library(ggplot2)
library(vcfR)
library(qqman)
library(pcadapt)


# First we use PCAdapt to read in our (uncompressed) VCF data from the class server 
vcf <- read.pcadapt("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf",
                    type="vcf")

# We also need the compressed version to subset the metadata file with
vcfR <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf.gz")

# NOTE: we ran into an issue in class and used the 2 files above that I had stored in the class directory
# If you wanted to use these on your own filtered VCF file, you can, you just need to provide paths to
# versions of the VCF file that are uncompressed (no .gz extension) in the read,pcadapt() line and to the
# compressed version (ends in .gz) in the read.vcfR() line
# Be careful to only save the uncompressed version if you make one outside your github repo, because these
# are big and would cause github to not accept your push

# Now, read in the meta data and subset it for the individuals in the vcfR object
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta2 <- meta[meta$id %in% colnames(vcfR@gt[,-1]),]

# Now run PCAdapt! Options annotated below...
pcadapt.pca <- pcadapt(vcf,
                       K=2,
                       method="componentwise",
                       min.maf=0.01,
                       LD.clumping = list(size=500, thr=0.2))

# K: number of PC axes you want to test for selection
# method: Selecting "componentwise" returns a separate test for selection on each PC axis; otherwise it returns a global selection test across all PC axes together
# min.maf: the minor allele frequency value to use when testing -- SNPs below this value won't get tested
# LD.clumping: removes loci in strong LD when initially fitting the PCs, then tests for selection using all loci; size is distance between adjacent SNPs in bp; thr = threshold value of the correlation coefficient (r^2, a measure of LD) between adjacent loci

summary(pcadapt.pca) # shows an overview of what's contained in the pcadapt.pca object

# A quick a dirty Manahattan plot, but without chromosomal info
plot(pcadapt.pca, K=2) 

# To make the Manhattan plot where we can actually see the outliers as a 
# function of chromosome and position, we need to get the SNP info out 
# of the "fix" region of the vcfR object

# We can get a quick view of that region by combining View(head()):
View(head(vcfR@fix)) 

# this shows us that the chromosome and position info are the frist 2 columns.  
# Let's grab that info out:
vcfR.fix <- as.data.frame(vcfR@fix[,1:2])

# We need to create a new variable that relabels each chromosome to a numeric value
# Here, we'll just keep the 8 major chromosomes and not the small scaffolds
chr.main <- unique(vcfR.fix$CHROM)[1:8]

# Make a dataframe by binding together the unique chromosome names and number them sequentially from 1 to 8
chrnum <- as.data.frame(cbind(chr.main, seq(1,8,1)))

# Now extract the p-values from the pcadapt object and then bind together with the chromosome info
# Note -- rows number be exactly the same here! 
# I.e., they must be the same # of loci to cbind() properly
Pval <- pcadapt.pca$pvalues
pcadapt.MHplot <- cbind(vcfR.fix, Pval)

# Now use a tidyr:left_join command to assign the chromosome numbers we made to all the loci
pcadapt.MHplot <- left_join(chrnum, pcadapt.MHplot, join_by(chr.main==CHROM))

# Create the "SNP" variable that qqman needs
pcadapt.MHplot <- pcadapt.MHplot %>%
  mutate(SNP=paste0(chr.main,"_",POS))

# Set variables as numeric that contain numbers
pcadapt.MHplot$V2 = as.numeric(pcadapt.MHplot$V2)
pcadapt.MHplot$POS = as.numeric(pcadapt.MHplot$POS)
pcadapt.MHplot$pPC1 = as.numeric(pcadapt.MHplot[,4])
pcadapt.MHplot$pPC2 = as.numeric(pcadapt.MHplot[,5])

# Drop loci with "NA" values -- these got filtered out by the min.maf option in the pcadapt() call
pcadapt.MHplot <- pcadapt.MHplot %>% drop_na(pPC1)

# Finally...the Manhattan plots!
# Let's make separate ones for selection outleirs on PC1 and then do one for PC2 next 
manhattan(pcadapt.MHplot,
          chr="V2",
          bp="POS",
          p="pPC1",
          col=c("blue4","orange3"),
          logP=T,
          ylab="-log10 p-value",
          genomewideline=F,
          main="PCAdapt genome scan for selection (PC1)")

manhattan(pcadapt.MHplot,
          chr="V2",
          bp="POS",
          p="pPC2",
          col=c("blue4","orange3"),
          logP=T,
          ylab="-log10 p-value",
          genomewideline=F,
          main="PCAdapt genome scan for selection (PC2)")

# There are more sophisticated ways to pull out the identities of the outliers
# but a quick way we did in class was to filter based on the quantile of the p-values,
# here, choosing a quantile of 0.001, corresponding to the smallest 0.1% of p-values as outliers

View(pcadapt.MHplot %>%
       filter(pPC1<quantile(pcadapt.MHplot$pPC1,0.001)))

# One could then use the NCBI genome browser for the Centaurea reference genome
# to look up the corresponding chromosome and position of the outliers and check
# what genes fall nearby. 
# https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCA_030169165.1

# For the outlier peak on chromosome 8 that we found
# it lands right in the middle of the MED14 gene, which is one of the protein 
# subunits of the Mediator complex that recruits RNApol-II to initiate transcription. 
# Studies in the model plant Arabidopsis thaliana have shown that MED14 is involved in 
# plant adaptation to cold and heat responses, as well as plant defense. 

# https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCA_030169165.1

# This offers an intriguing candidate gene involved in adaptation along the PC1 axis. Other
# approaches that use all the outlier loci en masse can be utilized to look for enrichment
# of certain functions, like gene ontology (GO) or biochemical pathways (KEGG).  We'll
# explore those in the next module on transcriptomics...






