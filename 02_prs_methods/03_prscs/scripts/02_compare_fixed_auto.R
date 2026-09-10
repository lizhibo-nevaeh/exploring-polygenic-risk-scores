library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

# GWAS summary statistics
gwas <- fread(
  file.path(
    work,
    "software/PRScs-master/test_data/sumstats_se.txt"
  )
)

# fixed phi = 0.01
fixed_file <- file.path(
  work,
  "results/official_chr22/eur_slurm_numpy2fixed_pst_eff_a1_b0.5_phi1e-02_chr22.txt"
)

fixed <- fread(fixed_file, header = FALSE)
setnames(
  fixed,
  c("CHR", "SNP", "BP", "A1_FIXED", "A2_FIXED", "BETA_FIXED")
)

# PRS-CS-auto
auto_file <- file.path(
  work,
  "results/official_chr22/eur_auto_numpy2fixed_v2_pst_eff_a1_b0.5_phiauto_chr22.txt"
)

auto <- fread(auto_file, header = FALSE)
setnames(
  auto,
  c("CHR", "SNP", "BP", "A1_AUTO", "A2_AUTO", "BETA_AUTO")
)

# Merge
dat <- merge(gwas, fixed, by = "SNP")
dat <- merge(dat, auto, by = "SNP")

# Align GWAS beta to the PRS-CS effect allele
dat[
  A1 == A1_FIXED & A2 == A2_FIXED,
  GWAS_BETA_ALIGNED := BETA
]

dat[
  A1 == A2_FIXED & A2 == A1_FIXED,
  GWAS_BETA_ALIGNED := -BETA
]

# Check allele consistency between fixed and auto runs
stopifnot(
  all(dat$A1_FIXED == dat$A1_AUTO),
  all(dat$A2_FIXED == dat$A2_AUTO)
)

dat[, SHRINK_FIXED :=
      abs(BETA_FIXED) / abs(GWAS_BETA_ALIGNED)]

dat[, SHRINK_AUTO :=
      abs(BETA_AUTO) / abs(GWAS_BETA_ALIGNED)]

cat("SNPs:", nrow(dat), "\n")

cat("\n===== Correlations =====\n")
cat(
  "GWAS vs fixed:",
  cor(dat$GWAS_BETA_ALIGNED, dat$BETA_FIXED),
  "\n"
)

cat(
  "GWAS vs auto :",
  cor(dat$GWAS_BETA_ALIGNED, dat$BETA_AUTO),
  "\n"
)

cat(
  "Fixed vs auto:",
  cor(dat$BETA_FIXED, dat$BETA_AUTO),
  "\n"
)

cat("\n===== Median absolute effects =====\n")
cat(
  "GWAS  :",
  median(abs(dat$GWAS_BETA_ALIGNED)),
  "\n"
)

cat(
  "Fixed :",
  median(abs(dat$BETA_FIXED)),
  "\n"
)

cat(
  "Auto  :",
  median(abs(dat$BETA_AUTO)),
  "\n"
)

cat("\n===== Median shrink ratio =====\n")
cat(
  "Fixed:",
  median(dat$SHRINK_FIXED, na.rm = TRUE),
  "\n"
)

cat(
  "Auto :",
  median(dat$SHRINK_AUTO, na.rm = TRUE),
  "\n"
)

top <- dat[
  order(-abs(GWAS_BETA_ALIGNED))
][
  1:15,
  .(
    SNP,
    GWAS_BETA_ALIGNED,
    BETA_FIXED,
    BETA_AUTO,
    SHRINK_FIXED,
    SHRINK_AUTO
  )
]

cat("\n===== Top 15 SNPs =====\n")
print(top)

fwrite(
  dat,
  file.path(work, "results/prscs_fixed_vs_auto_all.tsv"),
  sep = "\t"
)

fwrite(
  top,
  file.path(work, "results/prscs_fixed_vs_auto_top15.tsv"),
  sep = "\t"
)
