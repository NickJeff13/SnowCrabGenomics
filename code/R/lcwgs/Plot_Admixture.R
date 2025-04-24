library(dplyr)
library(ggplot2)

#read in .Q file  and combine results with .fam file
tbl=read.table("/mnt/sdb/SnowCrab_LCWGS/MAF_Filtered_Plink/snowcrab.maffiltered.2.Q")

#read .fam file

fam=read.table("data/admixture/exon.recode.fam")
fam$V1 <- gsub(".*i5.","",fam$V1)
fam$V1 <- gsub(".realigned.bam","",fam$V1)
fam$V2 <- gsub(".*i5.","",fam$V2)
fam$V2 <- gsub(".realigned.bam","",fam$V2)

#combine results and .fam file
results_q2=as.data.frame(cbind(fam[,1:2], tbl))
colnames(results_q2)=c("Pop", "ID", "q1","q2","q3") #"q4","q5","q6","q7","q8","q9","q10")

#Reorder individuals by Pop
results_q2[order(results_q2$Pop),]

#plot results quickly

ggplot(results_q2, aes(x=Individual))+
        

barplot(t(results_q2[order(results_q2$Pop),3:4]), 
        col=c("firebrick","blue","forestgreen", "gold", "turquoise", "grey","limegreen", "aliceblue", "orange","purple"),
        xlab="Individual #", ylab="Ancestry", border=NA, 
        names=results_q2[order(results_q2$Pop),1], las=2,space = 0, cex.names = 0.7)

