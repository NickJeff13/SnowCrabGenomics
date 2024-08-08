library(dplyr)
library(ggplot2)

#check CV
CV=read.table("CV_DU2_220K", header=F)
plot(CV$V4)

#read in .Q file  and combine results with .fam file
tbl=read.table("DU2_Lab_220K_nobadloc_maf005_noPR_thin200Kbp.3.Q")

#quick plot of results
barplot(t(as.matrix(tbl)), col=rainbow(3),
        xlab="Individual #", ylab="Ancestry", border=NA)

#read .fam file

fam=read.table("DU2_Lab_220K_nobadloc_maf005_noPR_thin200Kbp.fam")

#combine results and .fam file
results_q5=as.data.frame(cbind(fam[,1:2], tbl))
colnames(results_q5)=c("Pop", "ID", "q1", "q2", "q3")

#Reorder individuals by Pop
results_q5[order(results_q5$Pop),]

#plot results quickly
barplot(t(results_q5[order(results_q5$Pop),3:5]), 
        col=c("red", "dodgerblue4", "darkgoldenrod1"),
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

