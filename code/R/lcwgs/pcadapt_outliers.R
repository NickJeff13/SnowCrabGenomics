library(pcadapt)
library(qvalue)
library(dplyr)
library(ggplot2)
library(ggConvexHull)

setwd("/mnt/sdb/SnowCrab_LCWGS/")
all.snp <-read.pcadapt(input = "snowcrab.bed", type = "bed")

#choose k Principal components
y<-pcadapt(all.snp, K=3)
plot(y, option="manhattan")
plot(y, option="scores")

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
  plot(xx$loadings[,i], pch=19, cex=.3, ylab=paste0("Loadings PC",i))

# Maf and missingness filtered SNPs - 1082 inds, 8384961 SNPs
maf.filtered.snp <-read.pcadapt(input = "/mnt/sdb/SnowCrab_LCWGS/MAF_Filtered_Plink/snowcrab.maffiltered.bed", type = "bed")

maf.pcadapt <- pcadapt(maf.filtered.snp, K=6)

#read in inds file
maf.inds <- read_table("/mnt/sdb/SnowCrab_LCWGS/MAF_Filtered_Plink/snowcrab.maffiltered.fam", col_names = F)
head(maf.inds)
#extract sequencing batch info

seq.batch <- substr(maf.inds$X1, 1,7)

pca.with.depth <- cbind(b$scores, snp.depth)
pca.with.pops<-cbind(pca.with.depth, pops)

#various plots
plot(maf.pcadapt,option="screeplot")
plot(maf.pcadapt, option="scores", i=1, j=3)
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
pops <- c(rep("Mar C", 45),rep("Mar E",45),rep("Mar B",33),rep("Quebec",33), rep("Lilly Canyon",22), rep("CMA 3B",34),rep("CMA 6B",34), rep("CMA N5440", 21), rep("CMA 10B",6),"Baie des Chaleurs",rep("CMA 10B",7), "Baie des Chaleurs", rep("CMA 10B", 7), "Baie des Chaleurs", rep("CMA 10B", 7), "Baie des Chaleurs", rep("CMA 10B",7), "Baie des Chaleurs", rep("CMA 4", 7), "Baie des Chaleurs", rep("CMA 4", 7), "Baie des Chaleurs", rep("CMA 4", 7), "Baie des Chaleurs", rep("CMA 4", 7),"Baie des Chaleurs", rep("CMA 4", 6), rep("Baie des Chaleurs", 12), rep("CMA 5A", 33), rep("CMA 8A", 32), rep("Lilly Canyon", 9), rep("Baie des Chaleurs", 4), rep("CMA 3N200", 10), rep("Offshore GrandBanks", 33), rep("CMA 3D", 35), rep("CMA 10A", 34), rep("St. Marys Bay", 27), rep("Fortune Bay", 35), rep("Trinity Bay", 13), rep("West Cape Breton FEMALE", 8), rep("Bradelle Bank", 44), rep("CMA 3N200", 15), rep("Torngat", 3), "St. Marys Bay", rep("CMA N5440", 8), rep("West Cape Breton FEMALE", 4), rep("Trinity Bay", 20), rep("Northeast NS", 25), rep("NAFO 4X", 32), rep("Northeast NS Outer", 35), rep("CMA 12G", 33), rep("CMA 12C", 35), rep("LaurentianChannel", 35), rep("West Cape Breton FEMALE", 54), rep("Mar D", 8), rep("West Cape Breton MALE", 65), rep("Mar D", 26))

a <-pcadapt(exo.snps, K=10)

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
#plot with ggplot to colour points by sequencing depth 
pca.with.depth <- cbind(b$scores, snp.depth)
pca.with.pops<-cbind(pca.with.depth, pops)

#colour by sequencing depth
ggplot()+
  geom_point(data=pca.with.depth, aes(x=1, y=2, colour = MEAN_DEPTH))+
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

ggsave(filename = "PCAdapt_ExonSNPs_PC2_PC3_NoFacet_withHulls_byRegion.png", plot=last_plot(), path = "figures/", width=10, height=8, units="in", dpi=300)


#Try outlier detections on the exon derived SNPs
qvals<- qvalue(maf.pcadapt$pvalues)$qvalues
padj <- p.adjust(maf.pcadapt$pvalues, method="bonferroni")
alpha=0.001

outliers <-which(padj < alpha)
length(outliers)

#get the PCs associated with each outlier
snp_pc <-get.pc(maf.pcadapt, outliers)
