library(tidyverse)
library(readr)
library(ggbiplot)
library(Polychrome)

# read in data - start with the exome derived SNPs

pca <- read_table("data/plink pca/exon.eigenvec", col_names = FALSE)
eigenval <- scan("data/plink pca/exon.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
# spp
pca$ind <- gsub(".*i5.","",pca$ind)
pca$ind <- gsub(".realigned.bam","",pca$ind)

spp <- rep(NA, length(pca$ind))
spp[grep("^MC", pca$ind)] <- "Maritime C"
spp[grep("^ME", pca$ind)] <- "Maritime E"
spp[grep("^MB", pca$ind)] <- "Maritime B"
spp[grep("^Q", pca$ind)] <- "Quebec"
spp[grep("^LC", pca$ind)] <- "Lilly Canyon"
spp[grep("^3B", pca$ind)] <- "CMA 3B"
spp[grep("^6B", pca$ind)] <- "CMA 6B"
spp[grep("^TB", pca$ind)] <- "TrinityBay"
spp[grep("^10B", pca$ind)] <- "CMA 10B"
spp[grep("^4_", pca$ind)] <- "CMA 4"
spp[grep("^BdC", pca$ind)] <- "Chaleurs"
spp[grep("^5A", pca$ind)] <- "CMA 5A"
spp[grep("^8A", pca$ind)] <- "CMA 8A"
spp[grep("^CPS", pca$ind)] <- "NAFO 3L"
spp[grep("^3D", pca$ind)] <- "CMA 3D"
spp[grep("^10A", pca$ind)] <- "CMA 10A"
spp[grep("^SMB", pca$ind)] <- "StMarysBay"
spp[grep("^FB", pca$ind)] <- "FortuneBay"
spp[grep("^WCB_F", pca$ind)] <- "WCB_Female"
spp[grep("^WCB_M", pca$ind)] <- "WCB_Male"
spp[grep("^3N", pca$ind)] <- "CMA 3N200"
spp[grep("^NENS", pca$ind)] <- "NENS"
spp[grep("^4X", pca$ind)] <- "NAFO 4X"
spp[grep("^6222", pca$ind)] <- "CMA 12G"
spp[grep("^2420", pca$ind)] <- "CMA 12C"
spp[grep("^AW", pca$ind)] <- "LaurentianChannel"
spp[grep("^D0", pca$ind)] <- "Maritime D"
spp[grep("^BB", pca$ind)] <- "Bradelle"
spp[grep("^T[[:digit:]]", pca$ind)] <- "CMA N5440"




# location
loc <- rep(NA, length(pca$ind))
loc[grep("^MC", pca$ind)] <- "Maritime C"
loc[grep("^ME", pca$ind)] <- "Maritime E"
loc[grep("^MB", pca$ind)] <- "Maritime B"
loc[grep("^Q", pca$ind)] <- "Quebec"
loc[grep("^LC", pca$ind)] <- "Lily Canyon"
loc[grep("^3B", pca$ind)] <- "CMA 3B"
loc[grep("^6B", pca$ind)] <- "CMA 6B"
loc[grep("^TB", pca$ind)] <- "TrinityBay"
loc[grep("^10B", pca$ind)] <- "CMA 10B"
loc[grep("^4_", pca$ind)] <- "CMA 4"
loc[grep("^BdC", pca$ind)] <- "Chaleurs"
loc[grep("^5A", pca$ind)] <- "CMA 5A"
loc[grep("^8A", pca$ind)] <- "CMA 8A"
loc[grep("^CPS", pca$ind)] <- "NAFO 3L"
loc[grep("^3D", pca$ind)] <- "CMA 3D"
loc[grep("^10A", pca$ind)] <- "CMA 10A"
loc[grep("^SMB", pca$ind)] <- "StMarysBay"
loc[grep("^FB", pca$ind)] <- "FortuneBay"
loc[grep("^WCB_F", pca$ind)] <- "WCB_Female"
loc[grep("^WCB_M", pca$ind)] <- "WCB_Male"
loc[grep("^3N", pca$ind)] <- "CMA 3N200"
loc[grep("^NENS", pca$ind)] <- "NENS"
loc[grep("^4X", pca$ind)] <- "NAFO 4X"
loc[grep("^6222", pca$ind)] <- "CMA 12G"
loc[grep("^2420", pca$ind)] <- "CMA 12C"
loc[grep("^AW", pca$ind)] <- "LaurentianChannel"
loc[grep("^D0", pca$ind)] <- "Maritime D"
loc[grep("^BB", pca$ind)] <- "Bradelle"
loc[grep("^T[[:digit:]]", pca$ind)] <- "CMA N5440"

# combine - if you want to plot each in different colours
spp_loc <- paste0(spp, "_", loc)

# remake data.frame
pca <- as_tibble(data.frame(pca, spp, loc, spp_loc))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)


# plot pca
b <- ggplot(pca, aes(PC1, PC2)) + geom_point(aes(fill = spp), shape=21, col="black", size = 3) + stat_ellipse(aes(col=spp))
b <- b + scale_fill_manual(values = as.vector(glasbey.colors(n=29)))+facet_wrap(~loc) + scale_color_manual(values =as.vector(glasbey.colors(n=29)))
b <- b + coord_equal() + theme_bw()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

b + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

ggsave("Plink_Exon_PCA_facet_PC1PC2.png",plot = b, device = "png", path = "figures/", width = 16, height = 12, units = "in", dpi = 320)

#plot without facet 
b2 <- ggplot(pca, aes(PC1, PC2)) + geom_point(aes(fill = spp), shape=21, col="black", size = 3) +
  stat_ellipse(aes(col=spp))
b2 <- b2 + scale_fill_manual(values = as.vector(glasbey.colors(n=29))) + scale_color_manual(values =as.vector(glasbey.colors(n=29)))
b2 <- b2 + coord_equal() + theme_bw()
b2 + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

#b2 + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

ggsave("Plink_Exon_PCA_PC1PC2.png",plot = b2, device = "png", path = "figures/", width = 16, height = 12, units = "in", dpi = 320)

### Now do the full SNP dataset - 9 million SNPs!!
pca <- read_table("data/plink pca/snowcrab.eigenvec", col_names = FALSE)
eigenval <- scan("data/plink pca/snowcrab.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
# spp
pca$ind <- gsub(".*i5.","",pca$ind)
pca$ind <- gsub(".realigned.bam","",pca$ind)

spp <- rep(NA, length(pca$ind))
spp[grep("^MC", pca$ind)] <- "Maritime C"
spp[grep("^ME", pca$ind)] <- "Maritime E"
spp[grep("^MB", pca$ind)] <- "Maritime B"
spp[grep("^Q", pca$ind)] <- "Quebec"
spp[grep("^LC", pca$ind)] <- "Lilly Canyon"
spp[grep("^3B", pca$ind)] <- "CMA 3B"
spp[grep("^6B", pca$ind)] <- "CMA 6B"
spp[grep("^TB", pca$ind)] <- "TrinityBay"
spp[grep("^10B", pca$ind)] <- "CMA 10B"
spp[grep("^4_", pca$ind)] <- "CMA 4"
spp[grep("^BdC", pca$ind)] <- "Chaleurs"
spp[grep("^5A", pca$ind)] <- "CMA 5A"
spp[grep("^8A", pca$ind)] <- "CMA 8A"
spp[grep("^CPS", pca$ind)] <- "NAFO 3L"
spp[grep("^3D", pca$ind)] <- "CMA 3D"
spp[grep("^10A", pca$ind)] <- "CMA 10A"
spp[grep("^SMB", pca$ind)] <- "StMarysBay"
spp[grep("^FB", pca$ind)] <- "FortuneBay"
spp[grep("^WCB_F", pca$ind)] <- "WCB_Female"
spp[grep("^WCB_M", pca$ind)] <- "WCB_Male"
spp[grep("^3N", pca$ind)] <- "CMA 3N200"
spp[grep("^NENS", pca$ind)] <- "NENS"
spp[grep("^4X", pca$ind)] <- "NAFO 4X"
spp[grep("^6222", pca$ind)] <- "CMA 12G"
spp[grep("^2420", pca$ind)] <- "CMA 12C"
spp[grep("^AW", pca$ind)] <- "LaurentianChannel"
spp[grep("^D0", pca$ind)] <- "Maritime D"
spp[grep("^BB", pca$ind)] <- "Bradelle"
spp[grep("^T[[:digit:]]", pca$ind)] <- "CMA N5440"




# location
loc <- rep(NA, length(pca$ind))
loc[grep("^MC", pca$ind)] <- "Maritime C"
loc[grep("^ME", pca$ind)] <- "Maritime E"
loc[grep("^MB", pca$ind)] <- "Maritime B"
loc[grep("^Q", pca$ind)] <- "Quebec"
loc[grep("^LC", pca$ind)] <- "Lilly Canyon"
loc[grep("^3B", pca$ind)] <- "CMA 3B"
loc[grep("^6B", pca$ind)] <- "CMA 6B"
loc[grep("^TB", pca$ind)] <- "TrinityBay"
loc[grep("^10B", pca$ind)] <- "CMA 10B"
loc[grep("^4_", pca$ind)] <- "CMA 4"
loc[grep("^BdC", pca$ind)] <- "Chaleurs"
loc[grep("^5A", pca$ind)] <- "CMA 5A"
loc[grep("^8A", pca$ind)] <- "CMA 8A"
loc[grep("^CPS", pca$ind)] <- "NAFO 3L"
loc[grep("^3D", pca$ind)] <- "CMA 3D"
loc[grep("^10A", pca$ind)] <- "CMA 10A"
loc[grep("^SMB", pca$ind)] <- "StMarysBay"
loc[grep("^FB", pca$ind)] <- "FortuneBay"
loc[grep("^WCB_F", pca$ind)] <- "WCB_Female"
loc[grep("^WCB_M", pca$ind)] <- "WCB_Male"
loc[grep("^3N", pca$ind)] <- "CMA 3N200"
loc[grep("^NENS", pca$ind)] <- "NENS"
loc[grep("^4X", pca$ind)] <- "NAFO 4X"
loc[grep("^6222", pca$ind)] <- "CMA 12G"
loc[grep("^2420", pca$ind)] <- "CMA 12C"
loc[grep("^AW", pca$ind)] <- "LaurentianChannel"
loc[grep("^D0", pca$ind)] <- "Maritime D"
loc[grep("^BB", pca$ind)] <- "Bradelle"
loc[grep("^T[[:digit:]]", pca$ind)] <- "CMA N5440"

# combine - if you want to plot each in different colours
spp_loc <- paste0(spp, "_", loc)

# remake data.frame
pca <- as_tibble(data.frame(pca, spp, loc, spp_loc))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)


# plot pca
c <- ggplot(pca, aes(PC1, PC2)) + geom_point(aes(fill = spp), shape=21, col="black", size = 3)
c <- c + scale_fill_manual(values = as.vector(glasbey.colors(n=29)))+facet_wrap(~loc)
c <- c + coord_equal() + theme_bw()
c + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

c + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

ggsave("Plink_AllSNPs_PCA_facet_PC1PC2.png", plot = c, device = "png", 
       path = "figures/", width = 16, height = 12, units = "in", dpi = 320)

#non-faceted
e <- ggplot(pca, aes(PC1, PC2)) + geom_point(aes(fill = spp), shape=21, col="black", size = 3) #+ stat_ellipse(aes(col=spp))
e <- e + scale_fill_manual(values = as.vector(glasbey.colors(n=29)))
e <- e + coord_equal() + theme_bw()
e + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

e + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))
e + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

ggsave("Plink_AllSNPs_PCA_PC1PC2.png",plot = e, device = "png", path = "figures/", width = 16, height = 12, units = "in", dpi = 320)


##########
#Save Data
save.image("data/PlinkPCA.RData")
