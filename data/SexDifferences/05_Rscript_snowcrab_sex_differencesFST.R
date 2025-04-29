library(qqman)
library(ggplot2)
library(dplyr)

#Script for plotting results (barplot, pie chart) for top SNP that differentiated sexes 

#read snow crab data (plink .ped file for top SNP) and plot manhattan for sex differences
geno<-read.table("/filepath/snowcrab/SexFST/snowcrab_allinds_JACEEZ010007791_sexSNP.ped")

#Check data and combine alleles (V7/V8) to make genotype column
head(geno)
geno$genotype<-paste0(geno$V7, geno$V8)

#Load metadata for crabs
metad<-read.csv("/filepath/snowcrab/PCAdapt_results_100KSNPs_metadat.csv", header=T)
head(metad)

#combine genotype data and metadata
combined_dat<-merge(x=geno, by.x=2, y=metad, by.y=1)

#Subset data for sexed individuals only (males and females from WCB, BB)
fem<-combined_dat[which(combined_dat$Sex=="Female"),]
male<-combined_dat[which(combined_dat$Sex=="Male"),]
sexind<-rbind(fem, male)


#Plot sex and genotype data in barplot for individuals with known genotypes (exclude 00) and known sex:
ggplot(sexind[which(sexind$genotype!="00"),], aes(fill=Sex, x=genotype)) + 
  geom_bar(position='dodge', stat='count')+
  ggtitle("Genotypes of sexed invididuals at top SNP")+theme_classic()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15))#+facet_wrap(.~V1 ) #include facet to show pops seperately

#Plot same thing for all samples (exclude individuals with missing genotype, and include ONLY individuals with no sex info):
ggplot(combined_dat[which(combined_dat$genotype!="00" &
                            combined_dat$Sex==""),], aes(fill=Sex, x=genotype)) + scale_fill_manual(values="gray10")+
  geom_bar(position='dodge', stat='count')+
  ggtitle("Genotypes of unknown invididuals at top SNP")+theme_classic()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15), legend.position = "none")#+facet_wrap(.~V1 ) #include facet to show pops seperately 

#For all data get genotype frequencies for pops to make pie charts
#Summarize data - exclude individuals with missing genotype. 
prop_genotypedInd<-combined_dat[which(combined_dat$genotype!="00"),] %>%
  group_by(V1, genotype) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

#Pie chart of allele freq data
ggplot(data=prop_genotypedInd,
       aes(x="", y=freq, fill=genotype)) +
  scale_fill_manual(values=c("firebrick","orange","dodgerblue4"))+
  geom_bar(stat="identity", width=1) + 
  coord_polar(theta = "y") + 
  facet_wrap(.~V1, )  +theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())
