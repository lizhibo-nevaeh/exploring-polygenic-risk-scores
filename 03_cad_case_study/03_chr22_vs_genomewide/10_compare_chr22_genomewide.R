#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

chr22_file <- Sys.getenv("CHR22_PRS_FILE")
gw_file <- Sys.getenv("GENOMEWIDE_PRS_FILE")

outdir <- Sys.getenv("COMPARE_OUTDIR", unset = "results/final")
if (!nzchar(chr22_file) || !nzchar(gw_file)) stop("Set CHR22_PRS_FILE and GENOMEWIDE_PRS_FILE")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

chr22 <- fread(chr22_file)
gw <- fread(gw_file)

need22 <- c("IID", "PRS_Z")
needgw <- c("IID", "PRS_Z", "Population")

if (!all(need22 %in% names(chr22))) {
  stop("chr22 file missing columns: ",
       paste(setdiff(need22, names(chr22)), collapse=", "))
}
if (!all(needgw %in% names(gw))) {
  stop("genome-wide file missing columns: ",
       paste(setdiff(needgw, names(gw)), collapse=", "))
}

x <- merge(
  chr22[, .(IID, CHR22_PRS_Z = PRS_Z)],
  gw[, .(IID, GW_PRS_Z = PRS_Z, Population)],
  by = "IID"
)

if (nrow(x) == 0) stop("No overlapping IID between chr22 and genome-wide files.")

pearson <- cor(x$CHR22_PRS_Z, x$GW_PRS_Z, method="pearson")
spearman <- cor(x$CHR22_PRS_Z, x$GW_PRS_Z, method="spearman")
fit <- lm(GW_PRS_Z ~ CHR22_PRS_Z, data=x)
r2 <- summary(fit)$r.squared

stats <- data.table(
  Metric = c("Matched individuals", "Pearson_r", "Spearman_rho", "Linear_model_R2"),
  Value = c(nrow(x), pearson, spearman, r2)
)

anova_fit <- aov(GW_PRS_Z ~ Population, data=x)
anova_tab <- summary(anova_fit)[[1]]
anova_p <- anova_tab[["Pr(>F)"]][1]
kw <- kruskal.test(GW_PRS_Z ~ Population, data=x)

pop_stats <- data.table(
  Test = c("One-way ANOVA", "Kruskal-Wallis"),
  P_value = c(anova_p, kw$p.value)
)

fwrite(x, file.path(outdir, "CAD_chr22_vs_genomewide.tsv"), sep="\t")
fwrite(stats, file.path(outdir, "CAD_chr22_vs_genomewide_stats.tsv"), sep="\t")
fwrite(pop_stats, file.path(outdir, "CAD_genomewide_population_tests.tsv"), sep="\t")

png(
  file.path(outdir, "CAD_chr22_vs_genomewide.png"),
  width=1500, height=1200, res=170
)

plot(
  x$CHR22_PRS_Z, x$GW_PRS_Z,
  pch=16,
  xlab="Chromosome 22 CAD PRS (Z-score)",
  ylab="Genome-wide CAD PRS (Z-score)",
  main="chr22 pilot vs genome-wide CAD PRS"
)
abline(lm(GW_PRS_Z ~ CHR22_PRS_Z, data=x), lwd=2)

legend(
  "topleft",
  legend=c(
    sprintf("N = %d", nrow(x)),
    sprintf("Pearson r = %.3f", pearson),
    sprintf("R² = %.3f", r2)
  ),
  bty="n"
)

dev.off()

cat("===== chr22 vs genome-wide =====\n")
print(stats)

cat("\n===== population tests =====\n")
print(pop_stats)

cat("\nSaved:\n")
cat(file.path(outdir, "CAD_chr22_vs_genomewide.png"), "\n")
cat(file.path(outdir, "CAD_chr22_vs_genomewide_stats.tsv"), "\n")
cat(file.path(outdir, "CAD_genomewide_population_tests.tsv"), "\n")
