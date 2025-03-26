
# VCF To Treemix ----------------------------------------------------------

#Converting our VCFs to Treemix input format


# Load libraries ----------------------------------------------------------

library(dartR)
library(adegenet)


# Convert VCF to GenLight object ------------------------------------------

obj <- gl.read.vcf(vcffile = "/mnt/sdb/SnowCrab_LCWGS/allcrub.vcf.gz",
                   verbose = 3)

# now convert gl object to treemix
gl2treemix(obj, output="/mnt/sdb/SnowCrab_LCWGS/",
           verbose=2)