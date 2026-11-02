# Script to plot Admixture results from both PCANGSD and Admixture run on hard called genotypes


# Load libraries ----------------------------------------------------------
library(dplyr)
library(tidyverse)
library(ggplot2)


# LCWGS Admixture results -------------------------------------------------

#read in .Q file  and combine results with .fam file
tbl=read.table("data/admixture/snowcrab.maffiltered.3.Q")

#read .fam file

fam=read.table("data/admixture/snowcrab.maffiltered.fam")
fam$V1 <- gsub(".*i5.","",fam$V1)
fam$V1 <- gsub(".realigned.bam","",fam$V1)
fam$V2 <- gsub(".*i5.","",fam$V2)
fam$V2 <- gsub(".realigned.bam","",fam$V2)

#combine results and .fam file
results_q2=as.data.frame(cbind(fam[,1:2], tbl))
colnames(results_q2)=c("Pop", "ID", "q1","q2","q3") #"q4","q5","q6","q7","q8","q9","q10")

#Reorder individuals by Pop
results_q2[order(results_q2$Pop),]

#reshape to make it easier for ggplot

results_long <- results_q2 %>%
        pivot_longer(cols=c(q1, q2, q3), names_to = "membership", values_to = "value")

results_long$ID <- factor(results_long$ID, levels=unique(results_q2$ID))
#plot results quickly

ggplot(results_long, aes(x=ID, y=value, fill=membership))+
        geom_bar(stat="identity") +
        theme_bw()+
        theme(axis.text.x = element_text(angle=90, hjust=1))
        
# too many x-axis labels so let's thin them out a bit
breaks <- levels(results_long$ID)[seq(1, length(levels(results_long$ID)), by = 10)]
        
admix1 <- ggplot(results_long, aes(x = ID, y = value, fill = membership)) +
            geom_bar(stat = "identity", width=1) +
            #scale_x_discrete(breaks = breaks) +
            scale_x_discrete(labels = function(x) {
              ifelse(seq_along(x) %% 10 == 0, x, "")
            }
            )+
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            labs(x = "", y = "Proportion", fill = "Ancestry");admix1

ggsave(filename = "MAF_Filtered_Admixture_Q3.png", plot = last_plot(), device = "png", path = "figures/", width = 10, height=8, dpi = 300)

#old
# barplot(t(results_q2[order(results_q2$Pop),3:4]), 
#         col=c("firebrick","blue","forestgreen", "gold", "turquoise", "grey","limegreen", "aliceblue", "orange","purple"),
#         xlab="Individual #", ylab="Ancestry", border=NA, 
#         names=results_q2[order(results_q2$Pop),1], las=2,space = 0, cex.names = 0.7)
# 


# Plot Genotype Likelihoods Admixture -------------------------------------


tbl=read.table("data/pcangsd/snowcrub.big.admix.3.Q")

sites <- read.table("data/pcangsd/snowcrub.big.sites")
str(sites)

tbl$ind <- paste0("Ind_",1:1145)

tbl_long <- tbl %>%
  pivot_longer(cols=c(V1, V2, V3), names_to = "membership", values_to = "value")

breaks2 <- levels(tbl_long$ind)[seq(1, length(levels(tbl_long$ind)), by = 10)]


admix2 <- ggplot(tbl_long, aes(x=ind, y = value, fill = membership)) +
      geom_bar(stat = "identity", width=1) +
      #scale_x_discrete(breaks = breaks) +
  # scale_x_discrete(labels = function(x) {
  #   ifelse(seq_along(x) %% 10 ==0, x, "")
  # }
  # )+
      theme_minimal() +
      theme(axis.text.x = element_blank(), legend.position = "none") +
      labs(x = "", y = "Proportion", fill = "Ancestry");admix2


library(patchwork)

(admix2 + guides(fill="none"))/admix1 + 
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & 
  theme(legend.position = "bottom")

ggsave(filename = "figures/Admixture_LCWGS_and_PCANGSD_combined.png",
       plot = last_plot(),
       width = 12,
       height=8,
       dpi=300, device = "png")
