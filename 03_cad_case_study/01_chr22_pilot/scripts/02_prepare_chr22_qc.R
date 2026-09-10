library(data.table)

infile <- "data/gwas/CAD_chr22_GRCh37.tsv"

dat <- fread(infile)

# Allele harmonisation preparation
# Standardize alleles to uppercase
dat[, effect_allele := toupper(effect_allele)]
dat[, other_allele  := toupper(other_allele)]

# MAF = minor allele frequency
# Minor-allele frequency
dat[, MAF := pmin(
    effect_allele_frequency,
    1 - effect_allele_frequency
)]

valid_alleles <- c("A", "C", "G", "T")

# PRS-ready QC
# Retain:
# 1. A/C/G/T biallelic SNPs
# 2. Valid beta / SE / P values
# 3. MAF >= 1%
qc <- dat[
    effect_allele %in% valid_alleles &
    other_allele  %in% valid_alleles &
    effect_allele != other_allele &
    is.finite(beta) &
    is.finite(standard_error) &
    standard_error > 0 &
    is.finite(p_value) &
    p_value > 0 &
    p_value <= 1 &
    is.finite(MAF) &
    MAF >= 0.01
]

cat("Raw variants :", nrow(dat), "\n")
cat("PRS-QC SNPs  :", nrow(qc), "\n")

cat("\n===== N distribution =====\n")
print(summary(qc$n))

cat("\n===== Cases distribution =====\n")
print(summary(qc$cases))

# Controls = total N - cases
qc[, controls := n - cases]

# Effective sample size
# More appropriate than total N for imbalanced case-control data
qc[, N_eff := 4 / (1 / cases + 1 / controls)]

cat("\n===== Effective N distribution =====\n")
print(summary(qc$N_eff))

cat("\n===== MAF distribution =====\n")
print(summary(qc$MAF))

fwrite(
    qc,
    "data/gwas/CAD_chr22_QC_MAF01.tsv",
    sep = "\t"
)

cat("\nSaved:\n")
cat("data/gwas/CAD_chr22_QC_MAF01.tsv\n")
