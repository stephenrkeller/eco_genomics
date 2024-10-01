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


# Now, we will run Admixture analyses and create plots
# For Admixture, we're going to use the LEA R package.
# The function inside LEA is called 'snmf'

CentAdmix <- snmf("outputs/vcf_final.filtered.thinned.geno",
                  K=1:10,
                  entropy=T,
                  repetitions =3,
                  project="new")  # if you're adding to this analysis later, you could choose project="continue"

# We can compare evidence for different levels of K (or PCs) using 
# the cross-entropy from snmf and the screeplot from PCA:

par(mfrow=c(2,1)) # This sets up a multi-panel plot
plot(CentAdmix, col="blue4", main="SNMF") # This plots the Cross-Entropy score we can use for selecting models with K values that fit our data well
plot(CentPCA$eigenvalues[1:10], ylab="Eigenvalues", xlab="Number of PCs", col="blue4",main="PCA")
dev.off() # This turns off the multi-panel setting; need to do this, otherwise, all subsequent plots will be in 2 panels!

# Now, we set a value of "K" to investigate
myK=4

# Calculate the cross-entropy (=model fit; lower values are better) for all 
# reps, then determine which rep has the lowest score; we'll use that for plotting
CE = cross.entropy(CentAdmix, K=myK)
best = which.min(CE)

# Extract the ancestry coefficients (the Q scores)
myKQ = Q(CentAdmix, K=myK, run=best)

# and cbind to the metadata
myKQmeta = cbind(myKQ, meta2)

# set up a color panel to use
my.colors = c("blue4","gold","tomato","lightblue","olivedrab")

# sort the entire dataset by features of interest in the metadata prior to plotting
# Here, I first group by continent, then sort by region and pop within continents
myKQmeta  = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(region,pop, .by_group=TRUE)

# Lastly, make a ancestry plot, and save it as a pdf to my figures/ directory
pdf("figures/Amdixture_K4.pdf", width=10, height=5)
barplot(as.matrix(t(myKQmeta[,1:myK])),
        border=NA,
        space=0,
        col=my.colors[1:myK],
        xlab="Geographic regions",
        ylab="Ancestry proportions",
        main=paste0("Ancestry matrix K=",myK))
axis(1,
     at=1:length(myKQmeta$region),
     labels=myKQmeta$region,
     tick=F,
     cex.axis=0.5,
     las=3)
dev.off()

# You can create multiple plots like this by varying "myK=value" statement
# and then running the code below that again.  Everything should propagate from that myK value.


