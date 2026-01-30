library(ggplot2)
library(patchwork)


# Plotting Inbreeding coefficient from vcftools ---------------------------

#Exon derived SNPS
crab.het <- read.table("data/Exon_SNPs_InbreedingStats.csv", header = T, sep="\t") %>% glimpse()

crab.het <- crab.het %>%
  mutate(He=(N_SITES-E.HOM.)/N_SITES, Ho=(N_SITES-O.HOM.)/N_SITES, POP=gsub("Baie des Chal","Chaleurs",POP))

#plot Fis
p1 <- ggplot()+ 
  geom_boxplot(data = crab.het, aes(x=POP, y=F, fill=POP),colour="black",  alpha=0.7)+
  theme_bw()+
  #guides(fill=guide_legend(nrow=3,byrow=T))+
  theme(axis.text.x=element_text(angle=60, vjust=0.8, hjust=1), 
        text=element_text(size=18),
        legend.position = "none");p1

#Plot He/Ho
p2 <- ggplot()+ 
  geom_boxplot(data = crab.het, aes(x=POP, y=Ho, fill=POP),colour="black",  alpha=0.7)+
  geom_boxplot(data=crab.het, aes(x=POP, y=He), colour="black")+
  theme_bw()+
  xlab(label="Population")+
  #guides(fill=guide_legend(nrow=3,byrow=T))+
  theme(axis.text.x=element_text(angle=60, hjust=1),
        text=element_text(size=18),
        legend.position = "none");p2

#Full MAF-filtered dataset - 8,384,961 SNPs

all.snp.het <- read.table(file = "data/Full_SNPs_InbreedingStats.csv", header = T, sep = "\t") %>% glimpse()

all.snp.het <- all.snp.het %>% 
  mutate(He=(N_SITES-E.HOM.)/N_SITES, Ho=(N_SITES-O.HOM.)/N_SITES, POP=gsub("Baie des Chal","Chaleurs",POP))


p3 <- ggplot()+
  geom_boxplot(data=all.snp.het, aes(x=POP, y=F, fill=POP), colour="black", alpha=0.7)+
  theme_bw()+
  theme(text=element_text(size=18), 
        axis.text.x=element_blank(), 
        legend.position = "none",
        axis.title.x =element_blank());p3

p3/p1

ggsave("figures/CrabPops_FIS.png", plot = last_plot(), device = "png", width = 10, height=8,
       units = "in", dpi = 300)

p4 <- ggplot()+ 
  geom_boxplot(data = all.snp.het, aes(x=POP, y=Ho, fill=POP),colour="black",  alpha=0.7)+
  geom_boxplot(data=all.snp.het, aes(x=POP, y=He), colour="black")+
  theme_bw()+
  xlab(label ="")+
  #guides(fill=guide_legend(nrow=3,byrow=T))+
  theme(axis.text.x=element_blank(), 
        text=element_text(size=18),
        legend.position = "none");p4



# Tajima's D --------------------------------------------------------------

#Ran filtered all SNP dataset with Taj D in vcftools, 8,384,961 SNPs and 1145 inds

taj.d <- read.table("/mnt/sdb/SnowCrab_LCWGS/SummaryStats/allsnpsFiltered.Tajima.D", header = T) %>% glimpse()

taj.d <- taj.d %>% filter(!TajimaD == "NaN") %>% sample_n(10000)

taj.d$Position <- paste0(taj.d$CHROM,  sep="_", taj.d$BIN_START)


p5 <- ggplot()+
      geom_line(data=taj.d, aes(x=Position, y=TajimaD, group=1),colour="dodgerblue2")+
      labs(x= "SNP Position",
           y="Tajima's D")+
      theme(axis.text.x=element_blank(),
            text=element_text(size=18),
            axis.line = element_line(linewidth = 1));p5

ggsave("figures/AllSNPs_TajimasD.png", plot=p5, device="png", width=10, height=8, dpi=300)
      
p5/p4/p2


ggsave("figures/Combined_TajD_Heterozygosity.png",
       plot = last_plot(),
       device = "png",
       width=12,
       height=10)
