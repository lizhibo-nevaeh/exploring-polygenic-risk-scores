library(data.table)
x <- fread("results/CAD_chr22_1000G_EUR_PRS_standardized.tsv")

x[, Population := factor(Population,
 levels=c("CEU","FIN","GBR","IBS","TSI"))]

png("results/CAD_chr22_PRS_by_population.png",
    1400,1000,res=160)

boxplot(PRS_Z ~ Population, data=x,
        xlab="1000G EUR population",
        ylab="chr22 CAD PRS Z-score",
        main="Chromosome 22 CAD PRS by population")

abline(h=0,lty=2)

dev.off()
