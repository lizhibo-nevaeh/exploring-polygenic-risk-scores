library(bigsnpr)
library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
datadir <- file.path(work, "tmp-data")

bed <- file.path(datadir, "public-data3.bed")
backing <- file.path(work, "public_data3_bigsnp")
rds <- paste0(backing, ".rds")

# PLINK -> bigSNP
# Store backing files inside the project directory
if (!file.exists(rds)) {
  snp_readBed(
    bed,
    backingfile = backing
  )
}

obj <- snp_attach(rds)

G   <- obj$genotypes
map0 <- obj$map
fam <- obj$fam

sumstats <- fread(
  file.path(datadir, "public-data3-sumstats.txt")
)

cat("===== Target genotype =====\n")
cat("Individuals:", nrow(G), "\n")
cat("Variants:", ncol(G), "\n")

cat("\n===== External GWAS =====\n")
cat("Variants:", nrow(sumstats), "\n")
cat("GWAS N:", unique(sumstats$N), "\n")

cat("\n===== Phenotype =====\n")
print(summary(fam$affection))

# Continuous trait: n_eff = N
sumstats[, n_eff := N]

# Format the genotype map following the tutorial
map <- setNames(
  map0[-3],
  c("chr", "rsid", "pos", "a1", "a0")
)

cat("\n===== Variant harmonisation =====\n")

# Match by rsID for this tutorial dataset
df_beta <- snp_match(
  sumstats,
  map,
  join_by_pos = FALSE
)

cat("\nMatched variants:", nrow(df_beta), "\n")

saveRDS(
  df_beta,
  file.path(work, "matched_sumstats.rds")
)

fwrite(
  df_beta,
  file.path(work, "matched_sumstats.tsv"),
  sep = "\t"
)

cat("\n===== First matched variants =====\n")
print(head(df_beta))

cat("\nSaved matched_sumstats.rds / .tsv\n")

df_beta <- snp_match(
  sumstats,
  map,
  join_by_pos = FALSE
)

cat("\nMatched variants:", nrow(df_beta), "\n")

saveRDS(
  df_beta,
  file.path(work, "matched_sumstats.rds")
)

fwrite(
  df_beta,
  file.path(work, "matched_sumstats.tsv"),
  sep = "\t"
)

cat("\n===== First matched variants =====\n")
print(head(df_beta))

cat("\nSaved matched_sumstats.rds / matched_sumstats.tsv\n")
