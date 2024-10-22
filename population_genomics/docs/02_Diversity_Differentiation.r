# Estimating diversity and genetic differentiation in the filtered Centaurea data

library(vcfR)
library(tidyverse)
library(qqman)

# helps solve plotting issues:
X11.options(type="cairo")
options(bitmapType="cairo")

# read in our VCF file from out repo outputs/ directory

vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.0.5.vcf.gz")

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

# Sometimes its helpful to save out derived results like our genetic diversity and differentiation metrics
# so we can access them in the future without having to run through all the steps each time.
# Here, we'll save this as a csv file to our outputs directory
write.csv(vcf.div.MHplot, "~/projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion.csv",
          quote=F,
          row.names=F)

# OK, now on to looking at diversity *within* groups (Hs)
# Which columns in the vcf.div.MHplot house the Hs values for each region?
names(vcf.div.MHplot) #This prints them to screen, and we see the Hs values in columns 4 through 9

# Let's plot a histogram of those Hs values overall all the loci, 
# and overlay the results for each region in the same plot for comparison.
# We'll use tidy operations combined with the %>% operator to wrangel the data 
# into a form needed for plotting with ggplot.
vcf.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(position="identity", alpha=0.5, bins=50) +
  labs(title="Genome-wide expected heterozygosity (Hs)",fill="Regions",
       x="Gene diversity (Hs) within Regions", y="Counts of SNPs")

# We can save plots we live made with ggplot using the 'ggsave' command:
ggsave("Histogram_GenomDiversity_byRegion.pdf", 
       path="~/projects/eco_genomics/population_genomics/figures/")

# In addition to the histogram plots, it can be helpful to summarize the 
# results in terms of averages, variation (std dev) and sample size (N) 
# Tidy ehlps us accomplish this with the 'summarise' command, along with 
# preceding 'group_by' and 'filter' commands to get the comparisons we want.
# We can run the command without assigning the output to a new variable (as done in class)
# and it will just print results to the console, or we can assign it to a new variable 
# (here, named 'Hs_table')and then we can save that result to our repos using the write.csv command...
Hs_table <- vcf.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0 & value<0.25) %>%
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n())

write.csv(Hs_table, "population_genomics/outputs/Hs_table_noZeros.csv",
          quote=F,
          row.names=F)

# Note the above code was customized several times in class, based on our discussion of
# looking at Hs values across *all* loci, or just ones that had non-zero values (using the added filtering step)
# You're choice in any given analysis must be guided by what you want to learn, and it may be useful to look 
# at the data in several ways (i.e., with and without loci where Hs=0, or what proportion of loci are Hs=0 in each group).




