library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

w <- fread(
  file.path(work, "LDpred2_auto_weights.tsv")
)

# Absolute effect size
w[, ABS_GWAS := abs(GWAS_BETA)]
w[, ABS_LDPRED2 := abs(LDpred2_BETA)]

# Shrinkage ratio: LDpred2 weight / original GWAS weight
w[, SHRINK_RATIO :=
    fifelse(
      ABS_GWAS > 0,
      ABS_LDPRED2 / ABS_GWAS,
      NA_real_
    )
]

cat("Variants:", nrow(w), "\n")

cat("\n===== Overall comparison =====\n")
cat(
  "Correlation between GWAS beta and LDpred2 beta:",
  cor(w$GWAS_BETA, w$LDpred2_BETA),
  "\n"
)

cat(
  "Median |GWAS beta|:",
  median(w$ABS_GWAS),
  "\n"
)

cat(
  "Median |LDpred2 beta|:",
  median(w$ABS_LDPRED2),
  "\n"
)

cat(
  "Median shrink ratio:",
  median(w$SHRINK_RATIO, na.rm = TRUE),
  "\n"
)

cat("\n===== Top 15 SNPs by |GWAS beta| =====\n")

top <- w[
  order(-ABS_GWAS),
  .(
    SNP,
    GWAS_BETA,
    LDpred2_BETA,
    SHRINK_RATIO
  )
][1:15]

print(top)

fwrite(
  top,
  file.path(work, "top15_beta_comparison.tsv"),
  sep = "\t"
)
