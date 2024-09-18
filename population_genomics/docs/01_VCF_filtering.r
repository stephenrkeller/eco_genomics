library(vcfR)

# set working directory to the class data folder:
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

# Optional -- can use to show what's in the directory
list.files("variants/")

# use vcfR to read in your VCF data
vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")

# calling the object with no other options prints a quick summary
vcf

# Or can use the head() command to peek inside:
head(vcf)

# Initial visualizations can help show how diversity varies across the chromosome
# note, vcfR will bring in a single chromosome at a time for these viz plots

# this brings in the DNA sequence of the ref genome, in FASTA format:
dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

# this brings in the annotation file for the reference genome, showing where the genes are located:
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")

# For viz purposes, one can make a `chromR` object and then plot
chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)

plot(chr1) 

# The above plot shows histograms, but not how quality varies physically along the chromosome
# For that, we use the `chromoqc` plotting function. 
# Note the xlim() allows you to set the plot boundaries (in bp position) so you can zoom in/out on specific chromosome regions of interest

chromoqc(chr1, xlim=c(1e1, 1.1e8))

# OK, now that we've got some sense of the sequencing quality, let's look at the depth (DP)
# within the VCF file in the genotype (gt) data matrix where the actual data are:

DP <- extract.gt(vcf, element="DP", as.numeric=T) # this is a handy function to grab data from the VCF


DP[1:5,1:10] # we can print a subset of it to the screen to peek

quantile(DP) # we can also look at the frequency of DP -- lots of zeros, right? Let's set those to `NA`

DP[DP==0] <- NA 

quantile(DP, na.rm=T) # now e can look again and get a better sense of the distribution

# It's often handy to visualize this to get a big pitcture view
# we can do that with this heatmap.bp function, and can subset part of the DP matrix if it's too big...

heatmap.bp(DP[1:1000,], rlabels=F, clabels=F) #looks at just the 1st thousand SNPs

# OK, we've got some work to do with filtering
library(SNPfiltR) 

vcf.filt <- hard_filter(vcf, depth=3) # Could also explore more stringent DP values (say, 5, 10, ...)

# Filtering on max depth is also impt. Very high read count at a SNP can indicate mapping errors.
# Sometimes mapping errors are due to duplicated genes in the genome -- these should be mapped separately to
# their respective source genes, but sometimes if those gene duplicates (i.e. paralogs) aren't 
# divergent enough at the sequence level, then the reads will map to a single spot and hence read depth will 
# tend to double up.  that's the basis for a oft-used max depth filter of 2X the average DP/SNP.
vcf.filt <- max_depth(vcf.filt, maxdepth=60)

# The next filter requires having some metadata about our samples. Bring this in from the 
# course data share:
meta <- read.csv("metadata/meta4vcf.csv", header=T)

# Subset for 'id' and whatever column you may want to group samples by.
# Here, I choose the 'region' column (= col 4). 
# Note, you'll need to rename the grouping column to 'pop' since that's what the 
# Next function is expecting for the name of the 2nd column.
meta2 <- meta[,c(1,4)]

names(meta2) <- c("id","pop")

# Filters out individuals based on the amount of missing data at the individual-level
# How much should we tolerate? 75%? Look at the resulting plots to see how many
# individuals you lose when your % missing filter increases in stringency...
vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff=0.75)

# After all that filtering, some SNP sites are probably no longer polymorphic, 
# or they may have had >2 alleles from the start.  These 2 filters get rid of both:
# mac is minor allele count (i.e., how many copies of the `Alt` allele are observed?)

vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1) 

# Now we can do the same missingness filter step, but here
# we filter on SNP % missingness, (i.e., for a given SNP, how many of the 629 individuals are missing
# data at that SNP?)...what's a tolerable threshold for us?  A common threshold you'll see in the 
# literature is 50%.  In general, it's easier on a study to lose SNPs (we have tens of thousands) than individual samples.

vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff=0.5)


# that's the last step in our filtering.  We can now visualize the data matrix again
# by re-extracting the depth data for the SNPs and individuals that made it through out filtering
# and then plot with the heatmap again.
# Note we're going to call this `DP2` so we don't confuse it with the pre-filtered DP data above

DP2 <- extract.gt(vcf.filt.indSNPMiss,
                  element="DP",
                  as.numeric=T)

# Look at the quantile distribution:
quantile(DP2, na.rm=T)

# and plot as a heatmap...

heatmap.bp(DP2[1:5000,],
           rlabels=F, clabels=F)

#Once we're satisfied, we can save our new freshly filtered vcf file out
# to our repo for working with downstream so we don't have to repeat all these filtering steps each time
# Think carefully about *where* you're saving this file so you know how to get back to it... :)

write.vcf(vcf.filt.indSNPMiss, 
          "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")



