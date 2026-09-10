library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

auto <- readRDS(
  file.path(work, "ldpred2_auto_chains.rds")
)

# Chain QC metric used in the tutorial workflow
corr_range <- sapply(
  auto,
  function(x) diff(range(x$corr_est, na.rm = TRUE))
)

threshold <- 0.95 * quantile(
  corr_range,
  0.95,
  na.rm = TRUE
)

keep <- which(corr_range > threshold)

qc <- data.table(
  chain = seq_along(auto),
  corr_range = corr_range,
  keep = seq_along(auto) %in% keep
)

fwrite(
  qc,
  file.path(work, "ldpred2_auto_chain_qc.tsv"),
  sep = "\t"
)

cat("QC threshold:", threshold, "\n")
cat("Chains total:", length(auto), "\n")
cat("Chains kept :", length(keep), "\n")
cat("Chains kept IDs:\n")
print(keep)

cat("\nCorrelation range:\n")
print(qc)
