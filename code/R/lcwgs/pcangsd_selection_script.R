library(RcppCNPy)
cov <- as.matrix(read.table("data/pcangsd/snowcrub.big.cov")) # Reads estimated covariance matrix
sel <- as.matrix(read.table("data/pcangsd/snowcrub.big.selection")) # Reads PC based selection statistics

# Plot PCA plot
e <- eigen(cov)
plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd")

# Obtain p-values from PC-based selection scan
p <- pchisq(sel, 1, lower.tail=FALSE)
head(p)


qqchi<-function(x,...){
  lambda<-round(median(x)/qchisq(0.5,1),2)
  qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
  legend("topleft",paste("lambda=",lambda))
}

### read in seleciton statistics (chi2 distributed)

## make QQ plot to QC the test statistics
qqchi(sel)

# convert test statistic to p-value
pval<-1-pchisq(sel,1)

## read positions (hg38)
p<-read.table("data/pcangsd/snowcrub.big.sites")

names(p)<-c("chr","pos")

## make manhatten plot
plot(-log10(pval),xlab="Chromosomes",main="Manhattan plot")


## zoom into region
w<-range(which(pval<1e-7)) + c(-100,100)
keep<-w[1]:w[2]
plot(p$pos[keep],-log10(pval[keep]),col=p$chr[keep],xlab="HG38 Position chr2")

## see the position of the most significant SNP
p$pos[which.max(s)]