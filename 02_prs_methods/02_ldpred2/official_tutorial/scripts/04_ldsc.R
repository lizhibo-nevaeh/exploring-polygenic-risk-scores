library(bigsnpr)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

df_beta <- readRDS(
  file.path(work, "matched_sumstats_maf_qc.rds")
)

ld <- readRDS(
  file.path(work, "official_LD_scores.rds")
)

corr <- readRDS(
  file.path(work, "official_LD_sfbm.rds")
)

# Confirm that all three objects contain the same number of SNPs
stopifnot(
  nrow(df_beta) == length(ld),
  nrow(corr) == nrow(df_beta),
  ncol(corr) == nrow(df_beta)
)

cat("GWAS variants:", nrow(df_beta), "\n")
cat("LD scores:", length(ld), "\n")
cat("LD matrix:", nrow(corr), "x", ncol(corr), "\n")

# LD Score Regression
ldsc <- with(
  df_beta,
  snp_ldsc(
    ld,
    length(ld),
    chi2 = (beta / beta_se)^2,
    sample_size = n_eff,
    blocks = NULL
  )
)

cat("\n===== LD Score Regression =====\n")
print(ldsc)

cat("\nInitial h2 for LDpred2:", ldsc[["h2"]], "\n")

saveRDS(
  ldsc,
  file.path(work, "ldsc_result.rds")
)
