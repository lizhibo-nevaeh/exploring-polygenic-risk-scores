library(data.table)
chr <- as.integer(commandArgs(trailingOnly=TRUE)[1])
root <- Sys.getenv("GWPRS_ROOT"); hm3_file <- Sys.getenv("HM3_REF")
maf_min <- as.numeric(Sys.getenv("MAF_MIN")); maf_diff_max <- as.numeric(Sys.getenv("MAF_DIFF_MAX"))
raw <- fread(file.path(root,"data/gwas_split",sprintf("CAD_chr%d_raw.tsv",chr)))
req <- c("chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","cases","n")
stopifnot(all(req %in% names(raw)))
raw[, effect_allele:=toupper(effect_allele)]; raw[, other_allele:=toupper(other_allele)]
raw[, MAF:=pmin(effect_allele_frequency,1-effect_allele_frequency)]; raw[, controls:=n-cases]
raw[, N_eff:=4/(1/cases+1/controls)]; valid<-c("A","C","G","T")
qc <- raw[effect_allele %chin% valid & other_allele %chin% valid & effect_allele!=other_allele & is.finite(beta) & is.finite(standard_error) & standard_error>0 & is.finite(MAF) & MAF>=maf_min & cases>0 & controls>0]
ref <- fread(hm3_file)[CHR==chr]; setnames(ref,"MAF","REF_MAF"); ref[,`:=`(A1=toupper(A1),A2=toupper(A2))]
pos <- merge(qc,ref,by.x=c("chromosome","base_pair_location"),by.y=c("CHR","BP"),allow.cartesian=TRUE)
m <- copy(pos); m[,STATUS:=fifelse(effect_allele==A1 & other_allele==A2,"direct",fifelse(effect_allele==A2 & other_allele==A1,"swap","incompatible"))]
m <- m[STATUS!="incompatible"]; m[,BETA_ALIGNED:=fifelse(STATUS=="direct",beta,-beta)]; m[,MAF_DIFF:=abs(MAF-REF_MAF)]
m <- m[is.finite(MAF_DIFF) & MAF_DIFF<=maf_diff_max]; setorder(m,SNP,MAF_DIFF,standard_error); m <- m[!duplicated(SNP)]
out <- file.path(root,"data/harmonized"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
fwrite(m[,.(SNP,A1,A2,BETA=BETA_ALIGNED,SE=standard_error)],file.path(out,sprintf("CAD_chr%d_PRScs_sumstats.tsv",chr)),sep="\t")
fwrite(m[,.(chromosome,SNP,CM=0,base_pair_location,A1,A2)],file.path(out,sprintf("CAD_chr%d_final.bim",chr)),sep="\t",col.names=FALSE)
fwrite(unique(m[,.(chromosome,base_pair_location)]),file.path(out,sprintf("CAD_chr%d_positions.tsv",chr)),sep="\t",col.names=FALSE)
fwrite(m[,.(SNP)],file.path(out,sprintf("CAD_chr%d.snps",chr)),col.names=FALSE)
fwrite(m[,.(SNP,N_eff,MAF,REF_MAF,MAF_DIFF,STATUS)],file.path(out,sprintf("CAD_chr%d_meta.tsv",chr)),sep="\t")
s <- data.table(CHR=chr,RAW=nrow(raw),GWAS_QC=nrow(qc),POSITION_OVERLAP=nrow(pos),FINAL=nrow(m),MEDIAN_NEFF=median(m$N_eff,na.rm=TRUE))
fwrite(s,file.path(root,"results/qc",sprintf("harmonize_chr%d.tsv",chr)),sep="\t")
cat("CHR",chr,"FINAL SNPs:",nrow(m),"Median Neff:",round(median(m$N_eff)),"\n"); stopifnot(nrow(m)>0,!anyDuplicated(m$SNP))
