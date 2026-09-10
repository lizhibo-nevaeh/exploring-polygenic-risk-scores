library(data.table)
root<-Sys.getenv("GWPRS_ROOT"); chr<-as.integer(Sys.getenv("CHR"))
s<-fread(file.path(root,"data/harmonized",sprintf("CAD_chr%d_PRScs_sumstats.tsv",chr)))
b<-fread(file.path(root,"data/target",sprintf("chr%d/1000G_EUR_chr%d_CAD_HM3.bim",chr,chr)),col.names=c("CHR","SNP","CM","BP","B1","B2"))
x<-merge(s,b,by="SNP",all.x=TRUE); x[,STATUS:=fifelse(A1==B1 & A2==B2,"direct",fifelse(A1==B2 & A2==B1,"swap","incompatible"))]
fwrite(x[,.N,by=STATUS],file.path(root,"results/qc",sprintf("target_alleles_chr%d.tsv",chr)),sep="\t")
stopifnot(nrow(x)==nrow(s),!anyNA(x$B1),!anyDuplicated(b$SNP),sum(x$STATUS=="incompatible")==0); cat("chr",chr,"target allele QC: PASS\n")
