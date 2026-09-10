library(data.table)

x <- fread(
    "results/CAD_chr22_HM3_harmonized_full.tsv"
)

# Primary frequency QC
final <- x[
    is.finite(MAF_DIFF) &
    MAF_DIFF <= 0.10
]

excluded <- x[
    !is.finite(MAF_DIFF) |
    MAF_DIFF > 0.10
]

cat("Before frequency QC :", nrow(x), "\n")
cat("After frequency QC  :", nrow(final), "\n")
cat("Excluded            :", nrow(excluded), "\n")

cat("\n===== Final N =====\n")
print(summary(final$n))

cat("\n===== Final effective N =====\n")
print(summary(final$N_eff))

cat("\n===== Excluded SNPs =====\n")
print(
    excluded[
        order(-MAF_DIFF),
        .(
            SNP,
            base_pair_location,
            MAF,
            REF_MAF,
            MAF_DIFF,
            n,
            N_eff
        )
    ]
)

# PRS-CS input
out <- final[, .(
    SNP,
    A1,
    A2,
    BETA = BETA_ALIGNED,
    SE = standard_error
)]

fwrite(
    out,
    "data/gwas/CAD_chr22_PRScs_sumstats_final.tsv",
    sep = "\t"
)

fwrite(
    excluded,
    "results/CAD_chr22_frequency_outliers.tsv",
    sep = "\t"
)

cat("\nSaved final PRS-CS sumstats.\n")
