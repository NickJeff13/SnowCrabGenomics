library(pcadapt)
library(qvalue)
library(dplyr)
library(ggplot2)
library(ggConvexHull)
library(readr)
library(patchwork)
library(Polychrome)


setwd("/mnt/sdb/SnowCrab_LCWGS/")
all.snp <-read.pcadapt(input = "snowcrab.bed", type = "bed")

#choose k Principal components

x <- pcadapt(all.snp, K=10)
plot(x,option="screeplot")
plot(x, option="scores",i=1,j=2)
plot(x, option="scores", i=1, j=3)
plot(x, option="manhattan")

#reduce K=3 now
xx <- pcadapt(all.snp, K=3)
summary(xx)
plot(xx, option="qqplot")
plot(xx, option="stat.distribution")
plot(xx, option="manhattan")
plot(xx, option="scores")

par(mfrow=c(2,2))
for (i in 1:3)
  plot(xx$loadings[,i], pch=19, cex=.3, ylab=paste0("Loadings PC",i))+theme_bw()

# Maf and missingness filtered SNPs - 1082 inds, 8384961 SNPs
maf.filtered.snp <-read.pcadapt(input = "/mnt/sdb/SnowCrab_LCWGS/MAF_Filtered_Plink/snowcrab.maffiltered.bed", type = "bed")

maf.pcadapt <- pcadapt(maf.filtered.snp, K=5)

#read in inds file
maf.inds <- read_table("/mnt/sdb/SnowCrab_LCWGS/MAF_Filtered_Plink/snowcrab.maffiltered.fam", col_names = F)
head(maf.inds)
length(maf.inds$X1) #1082 inds

#extract population info
lst1 <- gregexpr('.', maf.inds$X1, fixed = TRUE)
pop.names1 <- substring(maf.inds$X1, sapply(lst1, `[`, 4) + 1, sapply(lst1, `[`, 5) - 1)
pop.names2 <- gsub("_[0-9]+", "", pop.names1)
#now lots of gsubbing 
pop.names2 <- gsub("MC.{1,2}","Mar C", pop.names2, ignore.case = F)
pop.names2 <- gsub("MC.{1,2}","Mar C", pop.names2, ignore.case = F)
pop.names2 <- gsub("MC.{1,2}","Mar C", pop.names2, ignore.case = F)
pop.names2 <- gsub("MC.{1,2}","Mar C", pop.names2, ignore.case = F)
pop.names2 <- gsub("MC.{1,2}","Mar C", pop.names2, ignore.case = F)
pop.names2 <- gsub("MC.{1,2}","Mar C", pop.names2, ignore.case = F)

#extract sequencing batch info

seq.batch <- substr(maf.inds$X1, 1,7)

#various plots
plot(maf.pcadapt,option="screeplot")
plot(maf.pcadapt, option="scores", i=1, j=2)
plot(maf.pcadapt, option="manhattan")

maf.pca.batch <- as.data.frame(cbind(maf.pcadapt$scores, seq.batch))

#colour by sequencing batch
seq.batch.plot<- ggplot(data=maf.pca.batch, aes(x=as.numeric(V1), y=as.numeric(V2), fill=seq.batch))+
                     geom_point(shape=21, size=3, colour="black")+
                     scale_fill_brewer(palette = "Set1")+
                      labs(fill="Sequencing Batch",
                           x="PCA 1",
                           y="PCA 2")+
                      theme_bw();seq.batch.plot

seq.batch.plot2<- ggplot(data=maf.pca.batch, aes(x=as.numeric(V1), y=as.numeric(V3), fill=seq.batch))+
  geom_point(shape=21, size=3, colour="black")+
  scale_fill_brewer(palette = "Set1")+
  labs(fill="Sequencing Batch",
       x="PCA 1",
       y="PCA 3")+
  theme_bw();seq.batch.plot2

seq.batch.plot / seq.batch.plot2

ggsave(filename = "Sequencing_Batch_Effect_PCAdapt.png", plot = last_plot(), device = "png", path = "~/Documents/GitHub/SnowCrabGenomics/figures/", width = 12, height = 10, dpi = 320)

# Exome SNPs --------------------------------------------------------------

exo.snps <-read.pcadapt(input = "exon.recode.bed", type = "bed")
#1072 individuals 134087 snps

#also read in the snp depth data as an effect
snp.depth <- read.table("exon.recodedep.idepth", header = T) %>% glimpse()
str(snp.depth) #1145 inds with 134087 SNPs
str(exo.snps) # 1072 inds with 134087 SNPs

sum(grepl(".BdC",inds))
inds<-gsub(x = snp.depth$INDV, pattern = ".realigned.bam", replacement = "")
pops <- c(rep("Mar C", 45),rep("Mar E",45),rep("Mar B",33),rep("Quebec",33), rep("Lilly Canyon",22), rep("CMA 3B",34),rep("CMA 6B",34), rep("CMA N5440", 21), rep("CMA 10B",6),"Chaleurs",rep("CMA 10B",7), "Chaleurs", rep("CMA 10B", 7), "Chaleurs", rep("CMA 10B", 7), "Chaleurs", rep("CMA 10B",7), "Chaleurs", rep("CMA 4", 7), "Chaleurs", rep("CMA 4", 7), "Chaleurs", rep("CMA 4", 7), "Chaleurs", rep("CMA 4", 7),"Chaleurs", rep("CMA 4", 6), rep("Chaleurs", 12), rep("CMA 5A", 33), rep("CMA 8A", 32), rep("Lilly Canyon", 9), rep("Chaleurs", 4), rep("CMA 3N200", 10), rep("NAFO 3L", 33), rep("CMA 3D", 35), rep("CMA 10A", 34), rep("St Marys Bay", 27), rep("Fortune Bay", 35), rep("Trinity Bay", 13), rep("West Cape Breton FEMALE", 8), rep("Bradelle Bank", 44), rep("CMA 3N200", 15), rep("CMA N5440", 3), "St Marys Bay", rep("CMA N5440", 8), rep("West Cape Breton FEMALE", 4), rep("Trinity Bay", 20), rep("NENSout", 25), rep("NAFO 4X", 32), rep("NENSin", 35), rep("CMA 12G", 33), rep("CMA 12C", 35), rep("Laurentian Chan", 35), rep("West Cape Breton FEMALE", 54), rep("Mar D", 8), rep("West Cape Breton MALE", 65), rep("Mar D", 26))

region <- c(rep("Scotian Shelf", 45),rep("Scotian Shelf", 45),rep("Scotian Shelf",33),rep("Quebec",33), rep("Newfoundland",22), rep("Newfoundland",34),rep("Newfoundland",34), rep("Labrador", 21), rep("Newfoundland",6),"GSL",rep("Newfoundland",7), "GSL", rep("Newfoundland", 7), "GSL", rep("Newfoundland", 7), "GSL", rep("Newfoundland",7), "GSL", rep("Newfoundland", 7), "GSL", rep("Newfoundland", 7), "GSL", rep("Newfoundland", 7), "GSL", rep("Newfoundland", 7),"GSL", rep("Newfoundland", 6), rep("GSL", 12), rep("Newfoundland", 33), rep("Newfoundland", 32), rep("Newfoundland", 9), rep("GSL", 4), rep("Newfoundland", 10), rep("Newfoundland", 33), rep("Newfoundland", 35), rep("Newfoundland", 34), rep("Newfoundland", 27), rep("Newfoundland", 35), rep("Newfoundland", 13), rep("GSL", 8), rep("GSL", 44), rep("Newfoundland", 15), rep("Labrador", 3), "Newfoundland", rep("Labrador", 8), rep("GSL", 4), rep("Newfoundland", 20), rep("Scotian Shelf", 25), rep("Scotian Shelf", 32), rep("Scotian Shelf", 35), rep("Newfoundland", 33), rep("Newfoundland", 35), rep("Newfoundland", 35), rep("GSL", 54), rep("Scotian Shelf", 8), rep("GSL", 65), rep("Scotian Shelf", 26))

a <-pcadapt(exo.snps, K=4)

#check K scree plot
plot(a, option= "screeplot")
plot(a, option= "scores")
plot(a, option="scores", i=2, j=3)
plot(a, option="scores", i=3, j=4)
plot(a, option="manhattan")

b=a 
summary(b)
plot(b, option="manhattan")
plot(b, option="qqplot")
plot(b, option="scores", pop=pops)
#plot with ggplot to colour points by sequencing depth 
pca.with.depth <- cbind(b$scores, snp.depth)
pca.with.pops<-cbind(pca.with.depth, pops)
pca.with.regions <- cbind(pca.with.pops, region)
colnames(pca.with.regions) <- c("PC1","PC2","PC3","PC4","IND","N_SITES","MEAN_DEPTH","pops", "region")

#colour by sequencing depth
ggplot()+
  geom_point(data=pca.with.pops, aes(x=PC3, y=PC4, colour = MEAN_DEPTH))+
  scale_colour_continuous(type="viridis")+
  theme_bw()

ggsave(filename = "ExonSNPs_pcadapt_PC1_PC3_bySeqDepth.png", plot = last_plot(), path = "~/Documents/GitHub/SnowCrabGenomics/figures/", width = 10, height=8, units = "in", dpi=300)


#colour by pop
ggplot(data=pca.with.pops, aes(x=PC1, y=PC2, colour = pops))+
  geom_point()+
  facet_wrap(vars(pops),scales="free")+
  theme_bw()

ggsave(filename = "PCAdapt_ExonSNPs_PC1_PC2_FacetByPop.png", plot = last_plot(), path = "~/Documents/GitHub/SnowCrabGenomics/figures/", width=10, height = 8, units = "in", dpi=300)

#colour by pop - no facet
p20 <- ggplot()+
    geom_point(data=pca.with.regions, aes(x=PC1, y=PC2, fill=region),
               color="black", size=3, shape=21)+
    scale_fill_manual(values = c("#c43b3b", "#80c43b", "#3bc4c4", "#7f3bc4","gold")) + 
    #scale_fill_manual(values =as.vector(glasbey.colors(n=32)))+
    #geom_convexhull(aes(fill=pops, color=pops), alpha=0.1)+
  guides(fill="none")+
    theme_bw();p20

ggsave(filename = "PCAdapt_ExonSNPs_PC1_PC2_NoFacet_withHulls_byRegion2.png",
       plot=p20, path = "~/Documents/GitHub/SnowCrabGenomics/figures/", 
       width=10, height=8, units="in", dpi=300)

p21 <- ggplot()+
  geom_point(data=pca.with.regions, aes(x=PC1, y=PC3, fill=region),
             color="black", size=3, shape=21)+
  scale_fill_manual(values = c("#c43b3b", "#80c43b", "#3bc4c4", "#7f3bc4","gold")) + 
  #scale_fill_manual(values =as.vector(glasbey.colors(n=32)))+
  #geom_convexhull(aes(fill=pops, color=pops), alpha=0.1)+
  theme_bw();p21

p20 + p21

ggsave(filename = "PCAdapt_ExonSNPs_PC123merged.png",
       plot=last_plot(), path = "~/Documents/GitHub/SnowCrabGenomics/figures/", 
       width=10, height=8, units="in", dpi=300)

## Make a map of PC1 and PC2 mean values 

    map.df <- data.frame(pca.with.pops$PC1, pca.with.pops$PC2, pca.with.pops$PC3,pca.with.pops$pops) %>%
      group_by(pca.with.pops.pops) %>%
      summarise(
        mean_PC1=mean(pca.with.pops.PC1),
        mean_PC2=mean(pca.with.pops.PC2),
        mean_PC3=mean(pca.with.pops.PC3)
      )


# Overlapping Poolseq-LCWGS SNPs ------------------------------------------
    
    ovlp.snps <-read.pcadapt(input = "/mnt/sdb/SnowCrab_LCWGS/snowcrab.poolfiltered.bed", type = "bed")
    #1094 individuals 1,834,608 snps
    
    #also read in the snp depth data as an effect
    snp.depth <- read.table("/mnt/sdb/SnowCrab_LCWGS/vcfs/out.idepth", header = T) %>% glimpse()

    
    ovlp.inds <- read.table("/mnt/sdb/SnowCrab_LCWGS/snowcrab.poolfiltered.fam", header = F) %>% glimpse()
    ovlp.snp.depth <- left_join(ovlp.inds, snp.depth, by=c("V1"="INDV"))
    
    ovlp.inds$V1 <- gsub(".*i5.","",ovlp.inds$V1)
    ovlp.inds$V1 <- gsub(".realigned.bam","",ovlp.inds$V1)
    
    pop <- rep(NA, length(ovlp.inds$V1))
    pop[grep("^MC", ovlp.inds$V1)] <- "Mar C"
    pop[grep("^ME", ovlp.inds$V1)] <- "Mar E"
    pop[grep("^MB", ovlp.inds$V1)] <- "Mar B"
    pop[grep("^Q", ovlp.inds$V1)] <- "Quebec"
    pop[grep("^LC", ovlp.inds$V1)] <- "Lilly Canyon"
    pop[grep("^3B", ovlp.inds$V1)] <- "CMA 3B"
    pop[grep("^6B", ovlp.inds$V1)] <- "CMA 6B"
    pop[grep("^TB", ovlp.inds$V1)] <- "TrinityBay"
    pop[grep("^10B", ovlp.inds$V1)] <- "CMA 10B"
    pop[grep("^4_", ovlp.inds$V1)] <- "CMA 4"
    pop[grep("^BdC", ovlp.inds$V1)] <- "Chaleurs"
    pop[grep("^5A", ovlp.inds$V1)] <- "CMA 5A"
    pop[grep("^8A", ovlp.inds$V1)] <- "CMA 8A"
    pop[grep("^CPS", ovlp.inds$V1)] <- "NAFO 3L"
    pop[grep("^3D", ovlp.inds$V1)] <- "CMA 3D"
    pop[grep("^10A", ovlp.inds$V1)] <- "CMA 10A"
    pop[grep("^SMB", ovlp.inds$V1)] <- "StMarysBay"
    pop[grep("^FB", ovlp.inds$V1)] <- "FortuneBay"
    pop[grep("^WCB_F", ovlp.inds$V1)] <- "WCB_Female"
    pop[grep("^WCB_M", ovlp.inds$V1)] <- "WCB_Male"
    pop[grep("^3N", ovlp.inds$V1)] <- "CMA 3N200"
    pop[grep("^NENS", ovlp.inds$V1)] <- "NENS"
    pop[grep("^4X", ovlp.inds$V1)] <- "NAFO 4X"
    pop[grep("^6222", ovlp.inds$V1)] <- "CMA 12G"
    pop[grep("^2420", ovlp.inds$V1)] <- "CMA 12C"
    pop[grep("^AW", ovlp.inds$V1)] <- "LaurentianChannel"
    pop[grep("^D0", ovlp.inds$V1)] <- "Mar D"
    pop[grep("^BB", ovlp.inds$V1)] <- "Bradelle"
    pop[grep("^T[[:digit:]]", ovlp.inds$V1)] <- "CMA N5440"
    
    a <-pcadapt(ovlp.snps, K=4)
    
    #check K scree plot
    plot(a, option= "screeplot")
    plot(a, option= "scores")
    plot(a, option="scores", i=2, j=3)
    plot(a, option="scores", i=1, j=2)
    plot(a, option="manhattan")
    
    b=a 
    summary(b)
    plot(b, option="manhattan")
    plot(b, option="qqplot")
    plot(b, option="scores", pop=pop)
    #plot with ggplot to colour points by sequencing depth 
    pca.with.depth <- cbind(b$scores, ovlp.snp.depth)
    pca.with.pops<-cbind(pca.with.depth, pop) %>% 
      dplyr::select(-c(V2,V3,V4,V5,V6))
    #pca.with.regions <- cbind(pca.with.pops, region)
    colnames(pca.with.pops) <- c("PC1","PC2","PC3","PC4","IND","N_SITES","MEAN_DEPTH","pops")
    
    #colour by sequencing depth
    ggplot()+
      geom_point(data=pca.with.pops, aes(x=PC1, y=PC3, colour = MEAN_DEPTH))+
      scale_colour_continuous(type="viridis")+
      theme_bw()
    
    ggsave(filename = "PoolFilteredSNPs_pcadapt_PC1_PC2_bySeqDepth.png", plot = last_plot(), path = "~/Documents/GitHub/SnowCrabGenomics/figures/", width = 10, height=8, units = "in", dpi=300)
    
    
    #colour by pop
    ggplot(data=pca.with.pops, aes(x=PC1, y=PC2, colour = pops))+
      geom_point()+
      facet_wrap(vars(pops),scales="free")+
      theme_bw()
    
    ggsave(filename = "PCAdapt_ExonSNPs_PC1_PC2_FacetByPop.png", plot = last_plot(), path = "~/Documents/GitHub/SnowCrabGenomics/figures/", width=10, height = 8, units = "in", dpi=300)
    
    #colour by pop - no facet
    p40 <- ggplot()+
      geom_point(data=pca.with.regions, aes(x=PC1, y=PC2, fill=pops),
                 color="black", size=3, shape=21)+
      #scale_fill_manual(values = c("#c43b3b", "#80c43b", "#3bc4c4", "#7f3bc4","gold")) + 
      scale_fill_manual(values =as.vector(glasbey.colors(n=30)))+
      #geom_convexhull(aes(fill=pops, color=pops), alpha=0.1)+
      #guides(fill="none")+
      theme_bw();p40
    
    ggsave(filename = "PCAdapt_PoolFilteredSNPs_PC1_PC3_NoFacet.png",
           plot=p40, path = "~/Documents/GitHub/SnowCrabGenomics/figures/", 
           width=10, height=8, units="in", dpi=300)
    
    
  (p20 + p21) / p40
    
    ggsave(filename = "PCAdapt_ExonSNPs_PC123merged_withPoolFilteredPlot.png",
           plot=last_plot(), path = "~/Documents/GitHub/SnowCrabGenomics/figures/", 
           width=10, height=8, units="in", dpi=300)
    
    
    
    
    
############### #Try outlier detections on the exon derived SNPs
qvals<- qvalue(maf.pcadapt$pvalues)$qvalues
padj <- p.adjust(maf.pcadapt$pvalues, method="bonferroni")
alpha=0.001

outliers <-which(padj < alpha)
length(outliers)

#get the PCs associated with each outlier
snp_pc <-get.pc(maf.pcadapt, outliers)

#load SNP bim file to subset

bim.snps <- read.table("/mnt/sdb/SnowCrab_LCWGS/MAF_Filtered_Plink/snowcrab.maffiltered.bim") %>% glimpse()
