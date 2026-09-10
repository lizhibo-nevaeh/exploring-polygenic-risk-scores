library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
script_dir <- dirname(normalizePath(file_arg))
work <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

score <- fread(
    file.path(
        work,
        "results/CAD_chr22_1000G_EUR_PRS.sscore"
    )
)

# PLINK2 may label the first column as #FID
setnames(
    score,
    old = names(score)[1],
    new = "FID"
)

score[, PRS := SCORE1_SUM]

# Standardized PRS
score[, PRS_Z := as.numeric(scale(PRS))]

cat("===== Overall PRS =====\n")
print(summary(score$PRS))

cat("\n===== Standardized PRS =====\n")
print(summary(score$PRS_Z))

cat("\nN:", nrow(score), "\n")
cat("Mean:", mean(score$PRS), "\n")
cat("SD:", sd(score$PRS), "\n")

cat("\n===== ALLELE_CT =====\n")
print(summary(score$ALLELE_CT))

# Old 1000G EUR population metadata
pop_file <- Sys.getenv("POP_PSAM")
if (!nzchar(pop_file)) stop("Set POP_PSAM to the 1000G EUR population metadata file")

pop <- fread(pop_file)

setnames(pop, old = names(pop)[1], new = "IID")

meta <- pop[, .(
    IID,
    SuperPop,
    Population
)]

out <- merge(
    score,
    meta,
    by = "IID",
    all.x = TRUE
)

cat("\n===== Population counts =====\n")
print(out[, .N, by = Population][order(Population)])

cat("\n===== PRS by population =====\n")
pop_summary <- out[
    ,
    .(
        N = .N,
        Mean_PRS = mean(PRS),
        SD_PRS = sd(PRS),
        Mean_Z = mean(PRS_Z),
        Median_Z = median(PRS_Z)
    ),
    by = Population
][order(Population)]

print(pop_summary)

# Percentile
out[, PRS_PERCENTILE :=
    frank(PRS, ties.method = "average") /
    .N * 100
]

setorder(out, -PRS)

cat("\n===== Highest 10 chr22 PRS =====\n")
print(
    out[
        1:10,
        .(
            IID,
            Population,
            PRS,
            PRS_Z,
            PRS_PERCENTILE
        )
    ]
)

cat("\n===== Lowest 10 chr22 PRS =====\n")
print(
    out[
        (.N-9):.N,
        .(
            IID,
            Population,
            PRS,
            PRS_Z,
            PRS_PERCENTILE
        )
    ]
)

fwrite(
    out,
    file.path(
        work,
        "results/CAD_chr22_1000G_EUR_PRS_standardized.tsv"
    ),
    sep = "\t"
)

fwrite(
    pop_summary,
    file.path(
        work,
        "results/CAD_chr22_PRS_population_summary.tsv"
    ),
    sep = "\t"
)

cat("\nSaved standardized PRS results.\n")
