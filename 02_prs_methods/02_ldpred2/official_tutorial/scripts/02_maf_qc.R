library(bigsnpr)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

obj <- snp_attach(
  file.path(work, "public_data3_bigsnp.rds")
)

G <- obj$genotypes

df_beta <- readRDS(
  file.path(work, "matched_sumstats.rds")
)

# All target individuals
ind_row <- rows_along(G)

# Use already harmonized SNPs only
maf <- snp_MAF(
  G,
  ind.row = ind_row,
  ind.col = df_beta$`_NUM_ID_`,
  ncores = 1
)

maf_thr <- 1 / sqrt(length(ind_row))

cat("Individuals:", length(ind_row), "\n")
cat("Matched SNPs before MAF QC:", nrow(df_beta), "\n")
cat("MAF threshold:", maf_thr, "\n")

df_beta_qc <- df_beta[maf > maf_thr, ]

cat("SNPs after MAF QC:", nrow(df_beta_qc), "\n")

saveRDS(
  df_beta_qc,
  file.path(work, "matched_sumstats_maf_qc.rds")
)
