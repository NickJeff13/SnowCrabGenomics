## Code adapted from Sarah's, figuring out the proportions of sex per population based on SNP genotypes


# Load libraries ----------------------------------------------------------

library(qqman)
library(ggplot2)
library(dplyr)
library(sf)
library(patchwork)
library(scatterpie)


#Script for plotting results (barplot, pie chart) for top SNP that differentiated sexes 

#read snow crab data (plink .ped file for top SNP) and plot manhattan for sex differences
geno<-read.table("data/SexDifferences/Results/snowcrab_allinds_JACEEZ010007791_sexSNP.ped")

#Check data and combine alleles (V7/V8) to make genotype column
head(geno)
geno$genotype <- paste0(geno$V7, geno$V8)

#Load metadata for crabs
metad<-read.csv("data/SexDifferences/files/PCAdapt_results_100KSNPs_metadat.csv", header=T)
head(metad)

#combine genotype data and metadata
combined_dat<-merge(x=geno, by.x=2, y=metad, by.y=1)

#Subset data for sexed individuals only (males and females from WCB, BB)
fem<-combined_dat[which(combined_dat$Sex=="Female"),]
male<-combined_dat[which(combined_dat$Sex=="Male"),]
sexind<-rbind(fem, male)


#Plot sex and genotype data in barplot for individuals with known genotypes (exclude 00) and known sex:
p1 <- ggplot(sexind[which(sexind$genotype!="00"),], aes(fill=Sex, x=genotype)) + 
  geom_bar(position='dodge', stat='count', color="black")+
  #ggtitle("Genotypes of sexed individuals at top SNP")+
  ylab(label = "Count")+
  scale_fill_manual(values=c("#D64550","dodgerblue4"))+
  theme_classic()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        text = element_text(size = 20));p1#+facet_wrap(.~V1 ) #include facet to show pops separately

ggsave("Sexed_Inds_genotypes.png", plot = p1, device = "png", path = "figures/", width = 10, height = 8, dpi = 300)

#Plot same thing for all samples (exclude individuals with missing genotype, and include ONLY individuals with no sex info):
p2 <- ggplot(combined_dat[which(combined_dat$genotype!="00" &
                            combined_dat$Sex==""),], aes(fill=genotype, x=genotype))+
  scale_fill_manual(values=c("#D64550","#D8D7BF","dodgerblue4"))+
  #scale_fill_brewer(palette = "Set2")+
  geom_bar(position='dodge', stat='count', color="black")+
  theme_classic()+
  ylab(label = "Count")+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15), legend.position = "none");p2#+facet_wrap(.~V1 ) #include facet to show pops separately

ggsave("Sexed_UnknownIndividuals.png", plot = p2, path = "figures/", device = "png", width = 10, height = 8, dpi=300)

#For all data get genotype frequencies for pops to make pie charts
#Summarize data - exclude individuals with missing genotype. 
prop_genotypedInd<-combined_dat[which(combined_dat$genotype!="00"),] %>%
  group_by(V1, genotype) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  as.data.frame()

# Replace site codes in prop_genotypedInd with matching codes from crab_coords

prop_genotypedInd$V1 <- gsub("10A", "CMA 10A", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("10B", "CMA 10B", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("2420", "CMA 12C", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("2N5440\\(Torngat\\)", "CMA N5440", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("3B", "CMA 3B", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("3D", "CMA 3D", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("4X", "NAFO 4X", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("5A\\(Bonavista\\)", "CMA 5A", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("6222", "CMA 12G", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("6B", "CMA 6B", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("8A", "CMA 8A", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("BaiedeChaleur", "Chaleurs", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("BradelleBank", "Bradelle Bank", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("FortuneBay", "Fortune Bay", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("LillyCanyon", "Lilly Canyon", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("ManagementArea4", "CMA 4", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("MaritimeB", "Mar B", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("MaritimeC", "Mar C", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("MaritimeD", "Mar D", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("MaritimeE", "Mar E", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("MatthewMariePierII/CPS", "NAFO 3L", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("QuebecNorthShore", "Quebec", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("R557_21", "Laurentian Chan", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("St.Mary'sBay", "St Marys Bay", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("TrinityBay", "Trinity Bay", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("WesternCapeBreton", "West Cape Breton", x = prop_genotypedInd$V1)
prop_genotypedInd$V1 <- gsub("3N200\\(3N\\)", "CMA 3N200", x = prop_genotypedInd$V1)


#check we still have all pops
length(unique(prop_genotypedInd$V1))

#Pie chart of allele freq data
p3<- ggplot(data=prop_genotypedInd,
       aes(x="", y=freq, fill=genotype)) +
  scale_fill_manual(values=c("#D64550","#D8D7BF","dodgerblue4"))+
  geom_bar(stat="identity", width=1,colour="black") + 
  labs("none")+
  coord_polar(theta = "y") + 
  facet_wrap(.~V1, )  + 
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x= element_blank(),
        axis.title.y = element_blank(),
        panel.grid  = element_blank(),
        strip.background = element_rect(fill="white"));p3

ggsave(filename = "SexGenotype_PiePlot_ByPop.png", plot = p3, path = "figures/", width = 10, height=8, dpi=300)

# Make a map with pie plots on it -----------------------------------------

#Load map data from Map script
load("data/RData/crabmapdata.RData")

combined_genos <- prop_genotypedInd %>%
  left_join(crab_coords, by=c("V1"="SampleSite"))

combined_genos <- combined_genos %>% 
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2]
  )

combined_genos <- combined_genos[,-8]

#for making the pie chart
plot_df <- combined_genos %>%
  st_drop_geometry() %>% 
  arrange(V1, genotype) %>%
  group_by(V1) %>%
  mutate(
    start = 2 * pi * cumsum(lag(freq, default = 0)),
    end = 2 * pi * cumsum(freq)
  ) %>%
  ungroup() %>%
  as.data.frame

pie_df <- plot_df %>%
  tidyr::pivot_wider(
    id_cols = c(V1, lon, lat),
    names_from = genotype,
    values_from = freq,
    values_fill = 0
  )

### set plot extent
plot_extent <- crab_coords%>%
  st_buffer(100*2000)%>%
  st_bbox()

p4 <- ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=bathy, fill=NA)+
  scatterpie::geom_scatterpie(
    aes(x = lon, y = lat, group = V1),
    data = pie_df,
    colour="black",
    cols = c("GG", "GT", "TT"),
    pie_scale = 1
  ) +
  coord_sf(expand=F, xlim=plot_extent[c(1,3)],ylim=plot_extent[c(2,4)])+
  theme_bw() +
  scale_fill_manual(values=c("#D64550","#D8D7BF","dodgerblue4"))+
  #scale_fill_brewer(palette = "Set2")+
  theme(text = element_text(size = 20),
        legend.position = "bottom")

p4 <- p4 + guides(fill=guide_legend(title="Genotype"));p4


ggsave(filename = "CrabSex_Map_byGenotype.png", plot = p4, device = "png", path = "figures/", 
      width = 10, height = 8, dpi = 300)

  ##Now make a figure for publication with patchwork and multiple panels

p4 + p1/p2 + plot_layout(widths=c(3,1)) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size=20))+
  plot_layout(widths=c(2,0.5))

ggsave(filename = "CrabSexGenotypes_Combined.png", 
       plot = last_plot(), path = 'figures/',
       device = "png", 
       dpi=300, width=15, height=10)
