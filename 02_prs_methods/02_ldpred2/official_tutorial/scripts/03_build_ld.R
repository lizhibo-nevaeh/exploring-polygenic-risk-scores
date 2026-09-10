library(bigsnpr)
library(bigsparser)
library(Matrix)
library(bigparallelr)
bigparallelr::set_blas_ncores(1)
options(default.nproc.blas = NULL)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

obj <- snp_attach(
  file.path(work, "public_data3_bigsnp.rds")
)

G <- obj$genotypes
POS2 <- obj$map$genetic.dist

df_beta <- readRDS(
  file.path(work, "matched_sumstats_maf_qc.rds")
)

cat("Variants for LD:", nrow(df_beta), "\n")

backing <- file.path(work, "official_LD_sfbm")

ld <- NULL
corr <- NULL

for (chr in 1:22) {

  # SNPs from the current chromosome
  ind.chr <- which(df_beta$chr == chr)

  if (length(ind.chr) == 0) next

  # Corresponding genotype-matrix columns
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]

  cat(
    "Chr", chr,
    "- SNPs:", length(ind.chr),
    "\n"
  )

  # 3 cM LD window
  corr0 <- snp_cor(
    G,
    ind.col   = ind.chr2,
    size      = 3 / 1000,
    infos.pos = POS2[ind.chr2],
    ncores    = 1
  )

  # LD score: sum of nearby r-squared values for each SNP
  ld_chr <- Matrix::colSums(corr0^2)

  if (is.null(corr)) {

    ld <- ld_chr

    corr <- as_SFBM(
      corr0,
      backing,
      compact = TRUE
    )

  } else {

    ld <- c(ld, ld_chr)

    corr$add_columns(
      corr0,
      nrow(corr)
    )
  }
}

corr$save()

saveRDS(
  ld,
  file.path(work, "official_LD_scores.rds")
)

cat("\n===== Genome-wide LD =====\n")
cat("Matrix dimensions:",
    nrow(corr), "x", ncol(corr), "\n")
cat("LD scores:", length(ld), "\n")
cat("SFBM size MB:",
    round(file.size(corr$sbk) / 1024^2, 2),
    "\n")

cat("Genome-wide LD matrix finished.\n")
