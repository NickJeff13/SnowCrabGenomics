library(dplyr)
library(ggplot2)

#check CV
CV=read.table("CV_DU2_220K", header=F)
plot(CV$V4)

#read in .Q file  and combine results with .fam file
tbl=read.table("data/admixture/exon.recode.10.Q")

#quick plot of results
barplot(t(as.matrix(tbl)), col=c("firebrick","blue","forestgreen", "gold", "turquoise", "grey","limegreen", "aliceblue", "orange","purple"),
        xlab="Individual #", ylab="Ancestry", border=NA)

#read .fam file

fam=read.table("data/admixture/exon.recode.fam")
fam$V1 <- gsub(".*i5.","",fam$V1)
fam$V1 <- gsub(".realigned.bam","",fam$V1)
fam$V2 <- gsub(".*i5.","",fam$V2)
fam$V2 <- gsub(".realigned.bam","",fam$V2)

#combine results and .fam file
results_q5=as.data.frame(cbind(fam[,1:2], tbl))
colnames(results_q5)=c("Pop", "ID", "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10")

#Reorder individuals by Pop
results_q5[order(results_q5$Pop),]

#plot results quickly
barplot(t(results_q5[order(results_q5$Pop),3:12]), 
        col=c("firebrick","blue","forestgreen", "gold", "turquoise", "grey","limegreen", "aliceblue", "orange","purple"),
        xlab="Individual #", ylab="Ancestry", border=NA, 
        names=results_q5[order(results_q5$Pop),1], las=2,space = 0, cex.names = 0.7)


#saves PDF of results

pdf(file="Admixture_Q3_DU2.pdf", width = 13.40, height=3.56)
barplot(t(results_q5[order(results_q5$Pop),3:5]), 
        col=c("red", "dodgerblue4", "darkgoldenrod1"),
        xlab="Individual #", ylab="Ancestry", border=NA, 
        names=results_q5[order(results_q5$Pop),1], las=2,space = 0, cex.names = 0.7)
dev.off()

#check to see if PDF is saved in working directory

