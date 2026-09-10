library(data.table)

sst <- fread(
    "data/gwas/CAD_chr22_PRScs_sumstats_final.tsv"
)

bim <- fread(
    "data/target/1000G_EUR_chr22_CAD_HM3.bim",
    col.names = c("CHR","SNP","CM","BP","BIM_A1","BIM_A2")
)

x <- merge(
    sst,
    bim,
    by = "SNP",
    all.x = TRUE
)

x[, STATUS :=
    fifelse(
        A1 == BIM_A1 & A2 == BIM_A2,
        "direct",
        fifelse(
            A1 == BIM_A2 & A2 == BIM_A1,
            "swap",
            "incompatible"
        )
    )
]

cat("Sumstats SNPs :", nrow(sst), "\n")
cat("Target SNPs   :", nrow(bim), "\n")
cat("Merged SNPs   :", sum(!is.na(x$BIM_A1)), "\n")

cat("\n===== Allele status =====\n")
print(x[, .N, by = STATUS])

cat("\nDuplicate sumstats IDs:",
    sum(duplicated(sst$SNP)), "\n")

cat("Duplicate target IDs:",
    sum(duplicated(bim$SNP)), "\n")

stopifnot(
    nrow(sst) == 16431,
    nrow(bim) == 16431,
    sum(is.na(x$BIM_A1)) == 0,
    sum(x$STATUS == "incompatible") == 0
)

cat("\nTARGET ALLELE QC: PASS\n")
