library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

chains <- readRDS(
  file.path(work, "ldpred2_auto_chains.rds")
)

diag <- rbindlist(
  lapply(seq_along(chains), function(i) {

    x <- chains[[i]]

    data.table(
      chain = i,

      p_final = x$p_est,
      p_min   = min(x$path_p_est, na.rm = TRUE),
      p_max   = max(x$path_p_est, na.rm = TRUE),

      h2_final = x$h2_est,
      h2_min   = min(x$path_h2_est, na.rm = TRUE),
      h2_max   = max(x$path_h2_est, na.rm = TRUE),

      alpha_final = x$alpha_est,
      alpha_min   = min(x$path_alpha_est, na.rm = TRUE),
      alpha_max   = max(x$path_alpha_est, na.rm = TRUE),

      beta_sd = sd(x$beta_est, na.rm = TRUE),
      beta_nonzero = sum(abs(x$beta_est) > 1e-12, na.rm = TRUE)
    )
  })
)

fwrite(
  diag,
  file.path(work, "ldpred2_auto_diagnostics.tsv"),
  sep = "\t"
)

print(diag)
