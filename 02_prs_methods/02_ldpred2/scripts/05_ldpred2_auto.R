library(bigsnpr)
library(bigsparser)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

# GWAS input prepared in the previous step
input <- readRDS(
  file.path(work, "ldpred2_input.rds")
)

df_beta <- input$df_beta
h2_init <- input$h2_init

# LD matrix stored as an SFBM
corr <- readRDS(
  file.path(work, "discovery_LD_sfbm.rds")
)

cat("Variants:", nrow(df_beta), "\n")
cat("Initial h2:", h2_init, "\n")
cat("LD matrix:", nrow(corr), "x", ncol(corr), "\n")

# Initialize multiple chains with different p values
p_init <- exp(
  seq(log(1e-4), log(0.2), length.out = 10)
)

set.seed(20260908)

multi_auto <- snp_ldpred2_auto(
  corr,
  df_beta,
  h2_init = h2_init,
  vec_p_init = p_init,
  burn_in = 500,
  num_iter = 500,
  allow_jump_sign = FALSE,
  shrink_corr = 0.95,
  ncores = 1
)

saveRDS(
  multi_auto,
  file.path(work, "ldpred2_auto_chains.rds")
)

# Summarize each chain
chain_summary <- data.frame(
  chain = seq_along(multi_auto),
  p_init = sapply(multi_auto, function(x) x$p_init),
  p_est = sapply(multi_auto, function(x) x$p_est),
  h2_est = sapply(multi_auto, function(x) x$h2_est),
  alpha_est = sapply(multi_auto, function(x) x$alpha_est),
  corr_range = sapply(
    multi_auto,
    function(x) diff(range(x$corr_est, na.rm = TRUE))
  ),
  beta_NA = sapply(
    multi_auto,
    function(x) sum(is.na(x$beta_est))
  )
)

write.table(
  chain_summary,
  file.path(work, "ldpred2_auto_chain_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\n===== LDpred2-auto chains =====\n")
print(chain_summary)

cat("\nLDpred2-auto finished.\n")
