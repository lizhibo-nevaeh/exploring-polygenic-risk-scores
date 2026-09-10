library(bigsnpr)
library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

obj <- snp_attach(file.path(work, "example_bigsnp_project.rds"))
G <- obj$genotypes
map <- obj$map

split <- fread(file.path(work, "sample_split.tsv"))
discovery <- sort(split[SET == "DISCOVERY", INDEX])

# Keep discovery SNPs with MAF >= 2%
maf <- snp_MAF(G, ind.row = discovery)
ind_snp <- which(is.finite(maf) & maf >= 0.02)

cat("Individuals in discovery:", length(discovery), "\n")
cat("SNPs before QC:", ncol(G), "\n")
cat("SNPs after MAF QC:", length(ind_snp), "\n")

# Principal component analysis
svd <- snp_autoSVD(
  G,
  infos.chr = map$chromosome,
  infos.pos = map$physical.pos,
  ind.row = discovery,
  ind.col = ind_snp,
  k = 10,
  ncores = 1
)

PC <- sweep(svd$u, 2, svd$d, "*")
colnames(PC) <- paste0("PC", 1:10)

fwrite(
  data.table(INDEX = discovery, PC),
  file.path(work, "discovery_PCs.tsv"),
  sep = "\t"
)

saveRDS(
  list(discovery = discovery, ind_snp = ind_snp, PC = PC),
  file.path(work, "pca_info.rds")
)

cat("PCA finished.\n")
