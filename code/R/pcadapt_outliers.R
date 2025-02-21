library(pcadapt)
library(qvalue)
library(dplyr)
library(ggplot2)
library(ggConvexHull)

setwd("/mnt/sdb/SnowCrab_LCWGS/")
all.snp <-read.pcadapt(input = "snowcrab.recode.bed", type = "bed")

#choose k Principal components
y<-pcadapt(all.snp, K=2)
plot(y, option="manhattan")
x <- pcadapt(all.snp, K=10)
plot(x,option="screeplot")
plot(x, option="scores")
plot(x, option="scores", i=1, j=3)
plot(x, option="manhattan")
#reduce K=3 now
xx <- pcadapt(all.snp, K=3)
summary(xx)
plot(xx, option="qqplot")
plot(xx, option="stat.distribution")

par(mfrow=c(2,2))
for (i in 1:3)
  plot(xx$loadings[,i], pch=19, cex=.3, ylab=paste0("Loadings PC",i))



#Try the imputed exome snps
exo.snps <-read.pcadapt(input = "exon.recode.bed", type = "bed")

#also read in the snp depth data as an effect
snp.depth <- read.table("exon.recodedep.idepth", header = T) %>% glimpse()
str(snp.depth) #1145 inds with 134087 SNPs
str(exo.snps) # 1072 inds with 134087 SNPs

sum(grepl(".BdC",inds))
inds<-gsub(x = snp.depth$INDV, pattern = ".realigned.bam", replacement = "")
pops <- c(rep("MC", 45),rep("ME",45),rep("MB",33),rep("Quebec",33), rep("Lily Canyon",22), rep("CFA 3B",34),rep("CFA 6B",34), rep("Torngat", 21), rep("CFA 10B",6),"Baie des Chaleurs",rep("CFA 10B",7), "Baie des Chaleurs", rep("CFA 10B", 7), "Baie des Chaleurs", rep("CFA 10B", 7), "Baie des Chaleurs", rep("CFA 10B",7), "Baie des Chaleurs", rep("CFA 4", 7), "Baie des Chaleurs", rep("CFA 4", 7), "Baie des Chaleurs", rep("CFA 4", 7), "Baie des Chaleurs", rep("CFA 4", 7),"Baie des Chaleurs", rep("CFA 4", 6), rep("Baie des Chaleurs", 12), rep("CFA 5A", 33), rep("CFA 8A", 32), rep("Lily Canyon", 9), rep("Baie des Chaleurs", 4), rep("3N", 10), rep("CPS", 33), rep("3D", 35), rep("CFA 10A", 34), rep("St. Marys Bay", 27), rep("Fortune Bay", 35), rep("Trinity Bay", 13), rep("West Cape Breton FEMALE", 8), rep("Bradelle Bank", 44), rep("3N", 15), rep("Torngat", 3), "St. Marys Bay", rep("Torngat", 8), rep("West Cape Breton FEMALE", 4), rep("Trinity Bay", 20), rep("Northeast NS", 25), rep("NAFO 4X", 32), rep("Northeast NS Outer", 35), rep("Labrador 6222", 33), rep("Labrador 2420", 35), rep("AW", 35), rep("West Cape Breton FEMALE", 54), rep("Dpop", 8), rep("West Cape Breton MALE", 65), rep("Dpop", 26))

a <-pcadapt(exo.snps, K=20)

#check K scree plot
plot(a, option= "screeplot")
plot(a, option= "scores")
plot(a, option="scores", i=2, j=3)
plot(a, option="scores", i=1, j=3)
plot(a, option="manhattan")

#re-run with K=3
b <-pcadapt(exo.snps, K=3)
summary(b)
plot(b, option="manhattan")
plot(b, option="qqplot")
plot(b, option="scores", pop=pops)
#plot with ggplot to colour points by sequenching depth 
pca.with.depth <- cbind(b$scores, snp.depth)
pca.with.pops<-cbind(pca.with.depth, pops)

#colour by sequencing depth
ggplot()+
  geom_point(data=pca.with.depth, aes(x=PC1, y=PC3, colour = DEPTH))+
  scale_colour_continuous(type="viridis")+
  theme_bw()

ggsave(filename = "ExonSNPs_PC1_PC3_bySeqDepth.png", plot = last_plot(), path = "~/Documents/GitHub/SnowCrabGenomics/figures/", width = 10, height=8, units = "in", dpi=300)


#colour by pop
ggplot(data=pca.with.pops, aes(x=PC1, y=PC3, colour = pops))+
  geom_point()+
  facet_wrap(vars(pops),scales="free")+
  theme_bw()

ggsave(filename = "PCAdapt_ExonSNPs_PC1_PC3_FacetByPop.png", plot = last_plot(), path = "~/Documents/GitHub/SnowCrabGenomics/figures/", width=10, height = 8, units = "in", dpi=300)

#colour by pop - no facet
ggplot(data=pca.with.pops, aes(x=PC1, y=PC2, colour = pops))+
  geom_point()+
  geom_convexhull(aes(fill=pops, color=pops), alpha=0.1)+
  theme_bw()

ggsave(filename = "PCAdapt_ExonSNPs_PC1_PC2_NoFacet_withHulls.png", plot=last_plot(), path = "~/Documents/GitHub/SnowCrabGenomics/figures/", width=10, height=8, units="in", dpi=300)


#Try outlier detections on the exon derived SNPs
qvals<- qvalue(b$pvalues)$qvalues
padj <- p.adjust(b$pvalues, method="bonferroni")
alpha=0.05

outliers <-which(padj < alpha)
length(outliers)

#get the PCs associated with each outlier
snp_pc <-get.pc(b, outliers)
