library(bigsnpr)
library(bigstatsr)
library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

obj <- snp_attach(file.path(work, "example_bigsnp_project.rds"))
info <- readRDS(file.path(work, "pca_info.rds"))

G   <- obj$genotypes
map <- obj$map
fam <- obj$fam

discovery <- info$discovery
ind_snp   <- info$ind_snp
PC        <- info$PC

# 0 = control, 1 = case
y <- as.integer(fam$affection == 2)

cat("Discovery N:", length(discovery), "\n")
cat("Controls:", sum(y[discovery] == 0), "\n")
cat("Cases:", sum(y[discovery] == 1), "\n")

# Fit logistic regression for each SNP
# Adjust for PC1-PC10
gwas <- big_univLogReg(
  G,
  y01.train   = y[discovery],
  ind.train   = discovery,
  ind.col     = ind_snp,
  covar.train = PC,
  ncores      = 1
)

# Calculate P values from Z scores
gwas$p.value <- predict(gwas, log10 = FALSE)

n_ctrl <- sum(y[discovery] == 0)
n_case <- sum(y[discovery] == 1)

# Effective sample size for the case-control GWAS
n_eff <- 4 / (1 / n_ctrl + 1 / n_case)

sumstats <- data.table(
  CHR   = map$chromosome[ind_snp],
  POS   = map$physical.pos[ind_snp],
  SNP   = map$marker.ID[ind_snp],
  A1    = map$allele1[ind_snp],
  A2    = map$allele2[ind_snp],
  BETA  = gwas$estim,
  SE    = gwas$std.err,
  P     = gwas$p.value,
  OR    = exp(gwas$estim),
  N     = length(discovery),
  N_EFF = n_eff
)

sumstats <- sumstats[
  is.finite(BETA) &
  is.finite(SE) &
  is.finite(P)
]

setorder(sumstats, P)

fwrite(
  sumstats,
  file.path(work, "discovery_GWAS_sumstats.tsv"),
  sep = "\t"
)

cat("\nEffective N:", round(n_eff, 2), "\n")
cat("Variants tested:", nrow(sumstats), "\n")

cat("\n===== Top 10 GWAS associations =====\n")
print(head(sumstats, 10))
