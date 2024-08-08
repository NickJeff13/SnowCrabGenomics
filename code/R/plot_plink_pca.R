library(tidyverse)
library(ggbiplot)

# read in data
pca <- read_table("exon.eigenvec", col_names = FALSE)
eigenval <- scan("exon.eigenval")
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
spp[grep("^LC", pca$ind)] <- "LC"
spp[grep("^3B", pca$ind)] <- "CFA3B"
spp[grep("^6B", pca$ind)] <- "CFA 6B"
spp[grep("^TB", pca$ind)] <- "TrinityBay"
spp[grep("^10B", pca$ind)] <- "CFA10B"
spp[grep("^4_", pca$ind)] <- "CFA 4"
spp[grep("^BdC", pca$ind)] <- "Chaleurs"
spp[grep("^5A", pca$ind)] <- "CFA5A"
spp[grep("^8A", pca$ind)] <- "CFA8A"
spp[grep("^CPS", pca$ind)] <- "CPS"
spp[grep("^3D", pca$ind)] <- "3D"
spp[grep("^10A", pca$ind)] <- "CFA10A"
spp[grep("^SMB", pca$ind)] <- "StMarysBay"
spp[grep("^FB", pca$ind)] <- "FortuneBay"
spp[grep("^WCB_F", pca$ind)] <- "WCB_Female"
spp[grep("^WCB_M", pca$ind)] <- "WCB_Male"
spp[grep("^3N", pca$ind)] <- "3N"
spp[grep("^NENS", pca$ind)] <- "NENS"
spp[grep("^4X", pca$ind)] <- "4X"
spp[grep("^6222", pca$ind)] <- "LAB6222"
spp[grep("^2420", pca$ind)] <- "LAB2420"
spp[grep("^AW", pca$ind)] <- "AW"
spp[grep("^D0", pca$ind)] <- "Dpop"
spp[grep("^BB", pca$ind)] <- "BB"
spp[grep("^T[[:digit:]]", pca$ind)] <- "TTTTT"




# location
loc <- rep(NA, length(pca$ind))
loc[grep("^MC", pca$ind)] <- "Maritime C"
loc[grep("^ME", pca$ind)] <- "Maritime E"
loc[grep("^MB", pca$ind)] <- "Maritime B"
loc[grep("^Q", pca$ind)] <- "Quebec"
loc[grep("^LC", pca$ind)] <- "LC"
loc[grep("^3B", pca$ind)] <- "CFA3B"
loc[grep("^6B", pca$ind)] <- "CFA 6B"
loc[grep("^TB", pca$ind)] <- "TrinityBay"
loc[grep("^10B", pca$ind)] <- "CFA10B"
loc[grep("^4_", pca$ind)] <- "CFA 4"
loc[grep("^BdC", pca$ind)] <- "Chaleurs"
loc[grep("^5A", pca$ind)] <- "CFA5A"
loc[grep("^8A", pca$ind)] <- "CFA8A"
loc[grep("^CPS", pca$ind)] <- "CPS"
loc[grep("^3D", pca$ind)] <- "3D"
loc[grep("^10A", pca$ind)] <- "CFA10A"
loc[grep("^SMB", pca$ind)] <- "StMarysBay"
loc[grep("^FB", pca$ind)] <- "FortuneBay"
loc[grep("^WCB_F", pca$ind)] <- "WCB_Female"
loc[grep("^WCB_M", pca$ind)] <- "WCB_Male"
loc[grep("^3N", pca$ind)] <- "3N"
loc[grep("^NENS", pca$ind)] <- "NENS"
loc[grep("^4X", pca$ind)] <- "4X"
loc[grep("^6222", pca$ind)] <- "LAB6222"
loc[grep("^2420", pca$ind)] <- "LAB2420"
loc[grep("^AW", pca$ind)] <- "AW"
loc[grep("^D0", pca$ind)] <- "Dpop"
loc[grep("^BB", pca$ind)] <- "BB"
loc[grep("^T[[:digit:]]", pca$ind)] <- "TTTTT"
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
b <- ggplot(pca, aes(PC2, PC3, col = spp)) + geom_point(size = 3)
b <- b + scale_color_manual(values = as.vector(glasbey(n=29)))+facet_wrap(~loc)
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))




##########
#Save Data
save.image("data/PlinkPCA.RData")