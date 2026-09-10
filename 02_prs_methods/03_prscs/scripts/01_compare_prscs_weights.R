library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

sst <- fread(
  file.path(
    work,
    "software/PRScs-master/test_data/sumstats_se.txt"
  )
)

pst <- fread(
  file.path(
    work,
    "results/official_chr22/eur_slurm_numpy2fixed_pst_eff_a1_b0.5_phi1e-02_chr22.txt"
  ),
  header = FALSE
)

setnames(
  pst,
  c("CHR", "SNP", "BP", "A1_PRSCS", "A2_PRSCS", "PRSCS_BETA")
)

dat <- merge(
  sst,
  pst,
  by = "SNP"
)

# Align GWAS beta to PRS-CS A1
dat[
  A1 == A1_PRSCS & A2 == A2_PRSCS,
  GWAS_BETA_ALIGNED := BETA
]

dat[
  A1 == A2_PRSCS & A2 == A1_PRSCS,
  GWAS_BETA_ALIGNED := -BETA
]

dat[, SHRINK_RATIO :=
      abs(PRSCS_BETA) / abs(GWAS_BETA_ALIGNED)]

cat("Matched SNPs:", nrow(dat), "\n")
cat(
  "Allele-aligned SNPs:",
  sum(!is.na(dat$GWAS_BETA_ALIGNED)),
  "\n"
)

cat(
  "Correlation:",
  cor(
    dat$GWAS_BETA_ALIGNED,
    dat$PRSCS_BETA,
    use = "complete.obs"
  ),
  "\n"
)

cat(
  "Median |GWAS beta|:",
  median(abs(dat$GWAS_BETA_ALIGNED), na.rm = TRUE),
  "\n"
)

cat(
  "Median |PRS-CS beta|:",
  median(abs(dat$PRSCS_BETA), na.rm = TRUE),
  "\n"
)

cat(
  "Median shrink ratio:",
  median(dat$SHRINK_RATIO, na.rm = TRUE),
  "\n"
)

top <- dat[
  !is.na(GWAS_BETA_ALIGNED)
][
  order(-abs(GWAS_BETA_ALIGNED))
][
  1:15,
  .(
    SNP,
    A1_PRSCS,
    A2_PRSCS,
    GWAS_BETA_ALIGNED,
    PRSCS_BETA,
    SHRINK_RATIO
  )
]

cat("\n===== Top 15 SNPs =====\n")
print(top)

fwrite(
  top,
  file.path(work, "results/prscs_top15_beta_comparison.tsv"),
  sep = "\t"
)
