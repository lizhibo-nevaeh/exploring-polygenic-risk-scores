library(data.table)
x <- fread("results/CAD_chr22_1000G_EUR_PRS_standardized.tsv")

png("results/CAD_chr22_PRS_distribution.png",
    1400,1000,res=160)

hist(x$PRS_Z, breaks=25, freq=FALSE,
     xlab="Standardized chr22 CAD PRS (Z-score)",
     main="Chromosome 22 CAD PRS pilot")

lines(density(x$PRS_Z), lwd=2)
abline(v=0,lty=2)

dev.off()
