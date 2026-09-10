library(bigsnpr)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

df_beta <- readRDS(
  file.path(work, "matched_sumstats_maf_qc.rds")
)

corr <- readRDS(
  file.path(work, "official_LD_sfbm.rds")
)

ldsc <- readRDS(
  file.path(work, "ldsc_result.rds")
)

h2_init <- ldsc[["h2"]]

# LDpred2 inputs used here: beta, beta_se, and n_eff
beta_df <- data.frame(
  beta    = df_beta$beta,
  beta_se = df_beta$beta_se,
  n_eff   = df_beta$n_eff
)

# Thirty chains with different initial p values
p_init <- exp(
  seq(
    log(1e-4),
    log(0.2),
    length.out = 30
  )
)

set.seed(20260909)

auto <- snp_ldpred2_auto(
  corr,
  beta_df,
  h2_init = h2_init,
  vec_p_init = p_init,
  burn_in = 500,
  num_iter = 500,
  allow_jump_sign = FALSE,
  shrink_corr = 0.95,
  ncores = 1
)

saveRDS(
  auto,
  file.path(work, "ldpred2_auto_chains.rds")
)

summary_tab <- data.frame(
  chain = seq_along(auto),
  p_init = sapply(auto, function(x) x$p_init),
  p_est  = sapply(auto, function(x) x$p_est),
  h2_est = sapply(auto, function(x) x$h2_est),
  alpha_est = sapply(auto, function(x) x$alpha_est),
  beta_NA = sapply(
    auto,
    function(x) sum(is.na(x$beta_est))
  )
)

write.table(
  summary_tab,
  file.path(work, "ldpred2_auto_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(summary_tab)
