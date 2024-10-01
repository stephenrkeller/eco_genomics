# load the usual libraries, plus one new one -- LEA
library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA) # This is the Landscape & Ecological Association package

options(bitmapType="cairo") # Helps resolve plotting issues when they crop up

setwd("~/projects/eco_genomics/population_genomics/") # set your working directory 

vcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz") # bring in your original filtered vcf file

# We need to thin the SNPs for LD (linkage disequilibrium) before we run
# PCA and Admixture analyses to satisfy the assumptions of independence among loci
# We'll use the function from the SNPfiltR package, thinning to keep only 1 SNP per 500 bp
vcf.thin <- distance_thin(vcf, min.distance=500) 

# Now read in meta data like usual and subset to keep only individuals also
# present in the filtered and thinned vcf file - -this is borrowed code from our 02_Diversity script!
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

meta2 <- meta[meta$id %in% colnames(vcf@gt[, -1]) , ]

dim(meta2)

# We now have a thinned vcf file and the metadata to accompany it. For the next steps,
# we need to save this thinned vcf to file, extract it, and convert it to a different 
# format (.geno) needed by the LEA program.  This need to convert between formats is a 
# hassle, but a commonplace in bioinformatics, so get used to it ;)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

# We now need to uncompress this vcf file, but when we do it'll be too big for our
# repos (github caps individual files at something ridiculous like 25 Mb).
# So, we'll "hide" the uncompressed vcf file in our home directory (~/) but outside 
# of our repo...don't forget where it lives!  And if you make a mistake and save it within your
# repo, just be sure to move it outside your repo before you do a git commit > push.

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

# Now perform the conversion from uncompressed vcf to .geno format:
geno <- vcf2geno(input.file="/gpfs1/home/s/r/srkeller/vcf_final.filtered.thinned.vcf",
                 output.file="outputs/vcf_final.filtered.thinned.geno")

# Now we're ready to do the PCA!
CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE)

# Note to future self, if you've already done the PCA previously, 
# you can load the results in without running it again like so:
CentPCA <- load.pcaProject("vcf_final.filtered.thinned.pcaProject")

show(CentPCA)

plot(CentPCA)

# plot(CentPCA$projections,
#      col=as.factor(meta2$region))
# legend("bottomright", 
#        legend=as.factor(unique(meta2$region)),
#        fill=as.factor(unique(meta2$region)))

ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +
       geom_point(alpha=0.5) +
       labs(title="Centaurea genetic PCA",x="PC1 (2.2%)",y="PC2 (1.1%)",color="Region",shape="Continent") 
#        xlim(-10,10) +
#        ylim(-10,10)

ggsave("figures/CentPCA_PC1vPC2.pdf", width=6, height=6, units ="in")





