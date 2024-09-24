# Estimating diversity and genetic differentiation in the filtered Centaurea data

library(vcfR)
library(tidyverse)
library(qqman)

# helps solve plotting issues:
X11.options(type="cairo")
options(bitmapType="cairo")

# read in our VCF file from out repio outputs/ directory

vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

# read in our metadata - -info on population of origin, what region the pops come from, what continent, etc.
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta) # vcf files has 595 samples
dim(meta) # meta has 629 inds

# We need to subset the 'meta' dataframe so that 
# it has just those sample id's that are also present in our filtered vcf file.
# we'll do that using the special %in% operator, which finds matches between 2 objects:
# Here, we're matching the 'id' column in 'meta' with the sample id's in the vcf file
# Note the vcf file stores the sample id's as column names in the @gt slot...minus the first column which
# is reserved for the "INFO" field.

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),] 

# Check the dimensions to be sure it worked:
dim(meta2)  # Looks good!

# calculate diversity stats using the genetic_diff fxn in vcfR
vcf.div <- genetic_diff(vcf, 
                        pops=as.factor(meta2$region),
                        method="nei")

# look at the format of the dataframe using the 'str' command:
str(vcf.div)

# For plotting routines below, we need to assign each chromosome a numeric label instead of a long text label
# If we know how many chromosomes are in our data, we can simply number them from 1 to n.
# The next few lines do this by first finding how many unique chromosomes are present in the vcf file, 
# and then taking just the 1st 8 values (the "main" chromosomes and not the smaller fragmented scaffolds)
chr.main <- unique(vcf.div$CHROM)[1:8]

# We then use the 'seq' function to number the chromosomes from 1 to 8
chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

# And finally merge these numerically named chromosomes with the full diversity results:
vcf.div.MHplot <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))

# Filter out Gst values <0  and create a new "SNP" column that concatenates the chromosome ID with the bp position
vcf.div.MHplot <- vcf.div.MHplot %>%
                        filter(Gst>0) %>%
                        mutate(SNP=paste0(chr.main,"_",POS))

# Make these 2 fields numbers instead of characters:
vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2)
vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)

# and finally....plot!!!
manhattan(vcf.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp=F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999))

# Lots of options can be explored here.  The "suggestiveline" is just that...suggestive that those regions
# of the genome above the line are experiencing a very high degree of population differentiation.  Why you might
# ask...could be something inhibiting gene flow into that portion of the genome (like genes contributing to 
# reproductive isolation between species) or perhaps these are genomic regions experiencing selection for local adaptation
# that causes regions to diverge from each other.  

