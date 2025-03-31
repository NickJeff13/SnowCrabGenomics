library(ggplot2)
library(patchwork)


# Plotting Inbreeding coefficient from vcftools ---------------------------

#Exon derived SNPS
crab.het <- read.table("data/Exon_SNPs_InbreedingStats.csv", header = T, sep="\t") %>% glimpse()

p1 <- ggplot()+ 
  geom_boxplot(data = crab.het, aes(x=POP, y=F, fill=POP),colour="black",  alpha=0.7)+
  theme_bw()+
  guides(fill=guide_legend(nrow=3,byrow=T))+
  theme(axis.text.x=element_text(angle=45, vjust=0.8, hjust=1), legend.position = "bottom");p1



# all 8,384,961 SNPs

all.snp.het <- read.table(file = "data/Full_SNPs_InbreedingStats.csv", header = T, sep = "\t") %>% glimpse()

p2 <- ggplot()+
  geom_boxplot(data=all.snp.het, aes(x=POP, y=F, fill=POP), colour="black", alpha=0.7)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, vjust=0.8, hjust=1), legend.position = "none",
        axis.title.x =element_blank());p2

p2/p1

ggsave("figures/CrabPops_Fis.png", plot = last_plot(), device = "png", width = 10, height=8,
       units = "in", dpi = 300)


# Tajima's D --------------------------------------------------------------

#Ran filtered all SNP dataset with Taj D in vcftools, 8,384,961 SNPs and 1145 inds

taj.d <- read.table("/mnt/sdb/SnowCrab_LCWGS/allsnpsFiltered.Tajima.D", header = T) %>% glimpse()

taj.d <- taj.d %>% filter(!TajimaD == "NaN") %>% sample_n(10000)

taj.d$Position <- paste0(taj.d$CHROM,  sep="_", taj.d$BIN_START)


p3 <- ggplot()+
      geom_line(data=taj.d, aes(x=Position, y=TajimaD, group=1))+
      labs(x= "SNP Position",
           y="Tajima's D")+
      theme_bw()+
      theme(axis.text.x=element_blank())

ggsave("figures/AllSNPs_TajimasD.png", plot=p3, device="png", width=10, height=8, dpi=300)
      
