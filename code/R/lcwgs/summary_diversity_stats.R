library(ggplot2)

crab.het <- read.csv("data/Exon_SNPs_InbreedingStats.csv", header = T) %>% glimpse()

p1 <- ggplot()+ 
  geom_boxplot(data = crab.het, aes(x=POP, y=F, fill=POP),colour="black",  alpha=0.7)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, vjust=0.8, hjust=1))

ggsave("figures/CrabPops_Fis.png", plot = p1, device = "png", width = 10, height=8,
       units = "in", dpi = 300)