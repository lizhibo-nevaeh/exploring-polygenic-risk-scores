library(data.table)

gwas_file <- "data/gwas/CAD_chr22_QC_MAF01.tsv"
ref_file <- Sys.getenv("PRSCS_HM3_REF")
if (!nzchar(ref_file)) stop("Set PRSCS_HM3_REF to snpinfo_1kg_hm3")

# 1. GWAS
gwas <- fread(gwas_file)

gwas[, effect_allele := toupper(effect_allele)]
gwas[, other_allele  := toupper(other_allele)]

cat("GWAS QC SNPs:", nrow(gwas), "\n")

# 2. PRS-CS HM3 reference
ref <- fread(ref_file)
ref <- ref[CHR == 22]

setnames(ref, "MAF", "REF_MAF")

ref[, A1 := toupper(A1)]
ref[, A2 := toupper(A2)]

cat("HM3 chr22 SNPs:", nrow(ref), "\n")

# 3. Match by chromosome + position
m <- merge(
    gwas,
    ref,
    by.x = c("chromosome", "base_pair_location"),
    by.y = c("CHR", "BP"),
    allow.cartesian = TRUE
)

cat("Position-overlap rows:", nrow(m), "\n")
cat("Position-overlap HM3 SNPs:", uniqueN(m$SNP), "\n")

# Palindromic SNP
m[, PALINDROMIC :=
    (effect_allele == "A" & other_allele == "T") |
    (effect_allele == "T" & other_allele == "A") |
    (effect_allele == "C" & other_allele == "G") |
    (effect_allele == "G" & other_allele == "C")
]

comp <- c(
    A = "T",
    T = "A",
    C = "G",
    G = "C"
)

m[, STATUS := NA_character_]
m[, BETA_ALIGNED := NA_real_]

# 4A. Direct match
m[
    !PALINDROMIC &
    effect_allele == A1 &
    other_allele == A2,
    `:=`(
        STATUS = "direct",
        BETA_ALIGNED = beta
    )
]

# 4B. Allele swap
m[
    !PALINDROMIC &
    effect_allele == A2 &
    other_allele == A1,
    `:=`(
        STATUS = "swap",
        BETA_ALIGNED = -beta
    )
]

# 4C. Strand complement
m[
    !PALINDROMIC &
    unname(comp[effect_allele]) == A1 &
    unname(comp[other_allele]) == A2,
    `:=`(
        STATUS = "strand",
        BETA_ALIGNED = beta
    )
]

# 4D. Strand complement + swap
m[
    !PALINDROMIC &
    unname(comp[effect_allele]) == A2 &
    unname(comp[other_allele]) == A1,
    `:=`(
        STATUS = "strand_swap",
        BETA_ALIGNED = -beta
    )
]

# Keep allele-compatible SNPs
matched <- m[!is.na(BETA_ALIGNED)]

# Require one unique compatible record per HM3 SNP
matched[, MATCH_N := .N, by = SNP]

unique_match <- matched[MATCH_N == 1]

# Frequency consistency check
unique_match[, MAF_DIFF := abs(MAF - REF_MAF)]

cat("\n===== Match summary =====\n")
cat("Matched rows:", nrow(matched), "\n")
cat("Unique matched HM3 SNPs:", nrow(unique_match), "\n")
cat(
    "Palindromic position-overlap rows:",
    sum(m$PALINDROMIC),
    "\n"
)

cat("\n===== Match status =====\n")
print(
    unique_match[, .N, by = STATUS][order(-N)]
)

cat("\n===== MAF difference =====\n")
print(summary(unique_match$MAF_DIFF))

cat("\n===== Matched total N =====\n")
print(summary(unique_match$n))

cat("\n===== Matched effective N =====\n")
print(summary(unique_match$N_eff))

# Detailed harmonisation result
fwrite(
    unique_match,
    "results/CAD_chr22_HM3_harmonized_full.tsv",
    sep = "\t"
)

# PRS-CS-ready sumstats
# A1/A2 follow the PRS-CS HM3 reference
# beta has been aligned to A1
prscs <- unique_match[, .(
    SNP,
    A1,
    A2,
    BETA = BETA_ALIGNED,
    SE = standard_error
)]

fwrite(
    prscs,
    "data/gwas/CAD_chr22_PRScs_sumstats.tsv",
    sep = "\t"
)

cat("\nSaved:\n")
cat("results/CAD_chr22_HM3_harmonized_full.tsv\n")
cat("data/gwas/CAD_chr22_PRScs_sumstats.tsv\n")
