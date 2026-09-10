library(bigsnpr)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

obj  <- snp_attach(file.path(work, "example_bigsnp_project.rds"))
info <- readRDS(file.path(work, "pca_info.rds"))

G   <- obj$genotypes
map <- obj$map

discovery <- info$discovery
ind_snp   <- info$ind_snp

cat("Individuals used for LD:", length(discovery), "\n")
cat("SNPs used for LD:", length(ind_snp), "\n")
cat("Chromosomes:", unique(map$chromosome[ind_snp]), "\n")

# Tutorial data: calculate correlations within 1000 kb of each SNP
corr <- snp_cor(
  G,
  ind.row   = discovery,
  ind.col   = ind_snp,
  size      = 1000,
  infos.pos = map$physical.pos[ind_snp],
  ncores    = 1
)

saveRDS(
  corr,
  file.path(work, "discovery_LD_matrix.rds")
)

cat("\nLD matrix dimensions:",
    nrow(corr), "x", ncol(corr), "\n")

cat("Non-zero entries:", length(corr@x), "\n")

cat("\n===== First 5 x 5 correlations =====\n")
print(round(as.matrix(corr[1:5, 1:5]), 3))

cat("\nLD matrix saved.\n")
