library(bigsnpr)
library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

obj <- snp_attach(
  file.path(work, "public_data3_bigsnp.rds")
)

G <- obj$genotypes
y <- obj$fam$affection

df_beta <- readRDS(
  file.path(work, "matched_sumstats_maf_qc.rds")
)

auto <- readRDS(
  file.path(work, "ldpred2_auto_chains.rds")
)

# Apply the same chain QC
corr_range <- sapply(
  auto,
  function(x) diff(range(x$corr_est, na.rm = TRUE))
)

threshold <- 0.95 * quantile(
  corr_range, 0.95, na.rm = TRUE
)

keep <- which(corr_range > threshold)

cat("Chains kept:", length(keep), "\n")

# Average chains that pass QC
beta_auto <- rowMeans(
  sapply(
    auto[keep],
    function(x) x$beta_est
  )
)

# Test split used in the tutorial workflow
set.seed(1)
ind_val  <- sample(rows_along(G), 350)
ind_test <- setdiff(rows_along(G), ind_val)

# Calculate PRS for test individuals
prs <- big_prodVec(
  G,
  beta_auto,
  ind.row = ind_test,
  ind.col = df_beta$`_NUM_ID_`
)

r  <- cor(prs, y[ind_test])
r2 <- r^2

cat("Test individuals:", length(ind_test), "\n")
cat("Pearson r:", r, "\n")
cat("R2:", r2, "\n")

# Save individual PRS
out <- data.table(
  INDEX = ind_test,
  IID = obj$fam$sample.ID[ind_test],
  PHENO = y[ind_test],
  PRS = prs,
  PRS_z = as.numeric(scale(prs))
)

fwrite(
  out,
  file.path(work, "LDpred2_auto_test_PRS.tsv"),
  sep = "\t"
)

# Save final SNP weights
weights <- data.table(
  SNP = df_beta$rsid,
  GWAS_BETA = df_beta$beta,
  LDpred2_BETA = beta_auto
)

fwrite(
  weights,
  file.path(work, "LDpred2_auto_weights.tsv"),
  sep = "\t"
)

cat("Results saved.\n")
