library(qqman)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggtext)

#Script for plotting FST results for sex comparison (Manhattan plots)

#read snow crab data and plot manhattan for sex differences - data is FST comparison between sexes (175 individuals from WCB, BB)
#Fst values for all 9M snps
all_fst<-data.table::fread("data/SexDifferences/snowcrab_sex_maf001_updateID_fst.fst")
#check file
nrow(all_fst)
head(all_fst)

#Update chromosome names to numbers (for plotting) - this may change them a bit*
all_fst$Chr_num<-gsub(all_fst$CHR,pattern="JACEEZ01", replacement = "")
all_fst$Chr_num<-gsub(all_fst$Chr_num,pattern=".1", replacement = "")
all_fst$Chr_num<-as.numeric(all_fst$Chr_num)

#Histogram of FST
hist(all_fst$FST)

#Density plot of fst values
ggplot(all_fst, aes(x=FST)) + 
  geom_density(col="dodgerblue", linewidth=1.25)+  # Overlay with transparent density plot
  theme_classic()

ggsave("figures/Sex_SNP_Fst_Density.png", plot = last_plot(), device = "png", 
       width = 10, height = 8, dpi =300)

#Subset data for plotting manhattan - don't need to plot all 9M snps, plot top SNPs and subset of low value SNPs
#subset top snps
high_fst<-all_fst[which(all_fst$FST>=0.05),]
#subset lower group of snps (note could have duplicates)
low_fst<-all_fst[which(all_fst$FST<0.2),]
subsample<-low_fst[sample(nrow(low_fst), 40000),]

#Use subset of data for plotting (9M SNP will be too large)
fst_for_plot<-rbind(high_fst,subsample)%>%
  as.data.frame()

#Manhattan plot for full dataset

data_cum <- fst_for_plot |>
  group_by(CHR) |>
  summarise(max_bp = max(POS)) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(CHR, bp_add)

gwas_data <- fst_for_plot |>
  inner_join(data_cum, by = "CHR") |>
  mutate(bp_cum = POS + bp_add)

manhattan(x = fst_for_plot, chr = "Chr_num", xlab="Scaffold", ylab="FST",
          bp = "POS",p = "FST",logp = F, col = c("black"))

axis_set <- gwas_data |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))


  manhplot <- ggplot(gwas_data, 
                     aes(x = bp_cum, 
                         y = FST)) +
    geom_point(alpha = 0.75) +
    scale_x_continuous(
      label = axis_set$CHR,
      breaks = axis_set$center
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.5,1)) +
    scale_color_manual(values = rep(
      c("#276FBF", "#183059"),
      unique(length(axis_set$CHR))
   )) +
    scale_size_continuous(range = c(0.5, 3)) +
    labs(
      x = NULL,
      y = "FST"
    ) +
   theme_bw() +
    theme(
     legend.position = "none",
     panel.grid.major.x = element_blank(),
     panel.grid.minor.x = element_blank(),
     axis.title.y = element_markdown(),
     axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    )

#Manhattan plot for Chr JACEEZ010007791.1 with highest Fst values
manhattan(x = all_fst[which(all_fst$CHR=="JACEEZ010007791.1"),],
          chr = "Chr_num", xlab="Chr JACEEZ010007791.1", ylab="FST",cex=1.5,
          bp = "POS",p = "FST",logp = F, col = c("black"), ylim=c(-0.05,1))

#subset FST values for chr JACEEZ010007791.1 
chr_sex<-all_fst[which(all_fst$CHR=="JACEEZ010007791.1"),]

#look at top SNPs
high_fst[which(high_fst$FST>0.9),]
region_sex<-chr_sex[which(chr_sex$FST>0.3),]
summary(region_sex$POS)

#Write table with top 100 SNPs 
write.table(all_fst[order(all_fst$FST,decreasing = T),][1:100],
            "/filepath/snowcrab/SexFST/top_FST_sex.txt",
            quote = F, row.names = F, col.names = T, sep="\t")



######Identifying HIGH FST regions with higher number of SNPs

nrow(high_fst[which(high_fst$FST>0.1),])


#Look for Chromosomes with high number of SNPs and focus on those
highdiff<-high_fst[which(high_fst$FST>0.1),]
high_fst_chr<-as.data.frame(table(highdiff$CHR))
high_fst_chr[order(high_fst_chr$Freq, decreasing = T),][1:10,]
more_than20<-high_fst_chr[which(high_fst_chr$Freq>20),]

#Chr JACEEZ010025865.1 is of interest
ggplot()+
  geom_bar(aes(x=more_than20$Var1,y=more_than20$Freq),
           stat="identity")+theme_classic()+
  theme(axis.text.x = element_text(angle=90))+
  ggtitle("Chr with >20 SNPs with high FST(>0.1)")+xlab("Chr")+ylab("# SNPs")


manhattan(x = all_fst[which(all_fst$CHR=="JACEEZ010025865.1"),],
          chr = "Chr_num", xlab="Chr JACEEZ010025865.1", ylab="FST",cex=1.5,
          bp = "POS",p = "FST",logp = F, col = c("black"), ylim=c(-0.05,0.5))

manhattan(x = all_fst[which(all_fst$CHR=="JACEEZ010025865.1"),],
          chr = "Chr_num", xlab="Chr JACEEZ010025865.1",
          ylab="FST",cex=1.5,
          bp = "POS",p = "FST",logp = F, col = c("black"), ylim=c(-0.05,0.5),
          xlim=c(280000,380000))
abline(v = 300000)
abline(v = 357000)

#extract SNPs in this region to make PCA
