library(data.table)

infile <- "data/gwas/CAD_chr22_GRCh37.tsv"

dat <- fread(
    infile,
    select = c(
        "p_value",
        "chromosome",
        "base_pair_location",
        "effect_allele",
        "other_allele",
        "effect_allele_frequency",
        "odds_ratio",
        "beta",
        "standard_error",
        "markername",
        "cases",
        "effective_cases",
        "n"
    )
)

# Standardize alleles to uppercase
dat[, effect_allele := toupper(effect_allele)]
dat[, other_allele  := toupper(other_allele)]

# MAF = minor allele frequency
# Convert EAF to minor-allele frequency as min(EAF, 1-EAF)
dat[, MAF := pmin(
    effect_allele_frequency,
    1 - effect_allele_frequency
)]

valid_allele <- c("A", "C", "G", "T")

is_snp <- (
    dat$effect_allele %in% valid_allele &
    dat$other_allele  %in% valid_allele &
    dat$effect_allele != dat$other_allele
)

valid_stats <- (
    is.finite(dat$beta) &
    is.finite(dat$standard_error) &
    dat$standard_error > 0 &
    is.finite(dat$p_value) &
    dat$p_value > 0 &
    dat$p_value <= 1
)

valid_freq <- (
    is.finite(dat$effect_allele_frequency) &
    dat$effect_allele_frequency > 0 &
    dat$effect_allele_frequency < 1
)

# Check consistency between beta and log(OR)
ok_or <- is.finite(dat$odds_ratio) & dat$odds_ratio > 0

beta_or_diff <- abs(
    dat$beta[ok_or] - log(dat$odds_ratio[ok_or])
)

summary_tab <- data.table(
    Metric = c(
        "Total_variants",
        "Valid_ACGT_biallelic_SNPs",
        "Valid_beta_SE_P",
        "Valid_EAF",
        "MAF_ge_0.01",
        "MAF_ge_0.05",
        "Duplicate_markername",
        "Duplicate_position",
        "Beta_logOR_match_1e-6"
    ),
    Value = c(
        nrow(dat),
        sum(is_snp),
        sum(valid_stats),
        sum(valid_freq),
        sum(dat$MAF >= 0.01, na.rm = TRUE),
        sum(dat$MAF >= 0.05, na.rm = TRUE),
        sum(duplicated(dat$markername)),
        sum(duplicated(dat$base_pair_location)),
        sum(beta_or_diff < 1e-6, na.rm = TRUE)
    )
)

cat("===== CAD chr22 GWAS QC =====\n")
print(summary_tab)

cat("\n===== Sample size =====\n")
cat("N range:\n")
print(range(dat$n, na.rm = TRUE))

cat("\nCases range:\n")
print(range(dat$cases, na.rm = TRUE))

cat("\n===== Alleles =====\n")
print(
    dat[, .N, by = .(effect_allele, other_allele)][
        order(-N)
    ]
)

cat("\n===== MAF summary =====\n")
print(summary(dat$MAF))

cat("\n===== P-value summary =====\n")
print(summary(dat$p_value))

fwrite(
    summary_tab,
    "results/chr22_gwas_qc_summary.tsv",
    sep = "\t"
)
