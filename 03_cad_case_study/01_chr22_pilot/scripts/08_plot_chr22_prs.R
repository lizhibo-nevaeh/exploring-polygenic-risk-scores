library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

dat <- fread(
    file.path(
        work,
        "results/CAD_chr22_1000G_EUR_PRS_standardized.tsv"
    )
)

figdir <- file.path(work, "results/figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Figure 1: overall distribution
# -----------------------------

png(
    file.path(figdir, "CAD_chr22_PRS_distribution.png"),
    width = 1800,
    height = 1300,
    res = 180
)

par(
    mar = c(5, 5, 4, 2),
    las = 1
)

hist(
    dat$PRS_Z,
    breaks = 25,
    freq = FALSE,
    border = "white",
    xlab = "Standardized chr22 CAD PRS (Z-score)",
    ylab = "Density",
    main = "Chromosome 22 CAD PRS pilot\n1000 Genomes European individuals"
)

lines(
    density(dat$PRS_Z),
    lwd = 2
)

abline(
    v = 0,
    lty = 2,
    lwd = 2
)

mtext(
    paste0(
        "N = ", nrow(dat),
        " | 16,431 SNPs | PRS-CS-auto"
    ),
    side = 3,
    line = 0.2,
    cex = 0.85
)

dev.off()


# -----------------------------
# Figure 2: population comparison
# -----------------------------

pop_order <- c("CEU", "FIN", "GBR", "IBS", "TSI")

dat[, Population := factor(
    Population,
    levels = pop_order
)]

png(
    file.path(figdir, "CAD_chr22_PRS_by_population.png"),
    width = 1800,
    height = 1300,
    res = 180
)

par(
    mar = c(5, 5, 4, 2),
    las = 1
)

boxplot(
    PRS_Z ~ Population,
    data = dat,
    xlab = "1000 Genomes European population",
    ylab = "Standardized chr22 CAD PRS (Z-score)",
    main = "Chromosome 22 CAD PRS across European populations",
    outline = FALSE
)

stripchart(
    PRS_Z ~ Population,
    data = dat,
    method = "jitter",
    vertical = TRUE,
    add = TRUE,
    pch = 16,
    cex = 0.55
)

abline(
    h = 0,
    lty = 2
)

mtext(
    "CEU 99 | FIN 99 | GBR 91 | IBS 107 | TSI 107",
    side = 3,
    line = 0.2,
    cex = 0.85
)

dev.off()

cat("Saved figures:\n")
cat(file.path(figdir, "CAD_chr22_PRS_distribution.png"), "\n")
cat(file.path(figdir, "CAD_chr22_PRS_by_population.png"), "\n")
