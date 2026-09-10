library(bigsnpr)
library(bigsparser)
library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

obj  <- snp_attach(file.path(work, "example_bigsnp_project.rds"))
info <- readRDS(file.path(work, "pca_info.rds"))

sumstats <- fread(file.path(work, "discovery_GWAS_sumstats.tsv"))
corr <- readRDS(file.path(work, "discovery_LD_matrix.rds"))

# SNP order in the LD matrix
snp_order <- obj$map$marker.ID[info$ind_snp]

# Reorder GWAS summary statistics to the same SNP order
m <- match(snp_order, sumstats$SNP)

stopifnot(!anyNA(m))

ss <- sumstats[m]

stopifnot(all(ss$SNP == snp_order))

# Core LDpred2 inputs
df_beta <- data.frame(
  beta    = ss$BETA,
  beta_se = ss$SE,
  n_eff   = ss$N_EFF
)

# Convert the LD matrix to an on-disk SFBM
# Store the backing file inside the project directory
backing <- file.path(work, "discovery_LD_sfbm")

corr_sfbm <- as_SFBM(
  corr,
  backingfile = backing,
  compact = TRUE
)

corr_sfbm$save()

# Estimate SNP heritability with LD score regression for LDpred2-auto initialization
ldsc <- snp_ldsc2(
  corr_sfbm,
  df_beta,
  blocks = NULL
)

h2_init <- as.numeric(ldsc["h2"])

cat("SNP order matched:", all(ss$SNP == snp_order), "\n")
cat("Variants:", nrow(df_beta), "\n")
cat("Initial h2 estimate:", h2_init, "\n")

saveRDS(
  list(
    snp_order = snp_order,
    df_beta = df_beta,
    h2_init = h2_init
  ),
  file.path(work, "ldpred2_input.rds")
)
