# Real-world CAD case study

## Why start with chromosome 22?

The first real-data analysis intentionally used only chromosome 22. The goal was to validate the technical chain:

real GWAS → QC → HapMap3 harmonization → PRS-CS-auto → posterior weights → target genotype scoring.

The chr22 pilot ended with **16,431 SNPs** scored in **503 1000G EUR individuals**.

## Why scale to the full genome?

CAD is highly polygenic. A single-chromosome score is useful for debugging but should not be treated as a substitute for a genome-wide PRS.

The final pipeline therefore applied the validated workflow across chromosomes 1–22.

## Final genome-wide results

Population summary of standardized PRS:

| Population | N | Mean Z | Median Z |
|---|---:|---:|---:|
| FIN | 99 | +0.467 | +0.418 |
| GBR | 91 | +0.153 | +0.090 |
| CEU | 99 | -0.144 | -0.146 |
| TSI | 107 | -0.176 | -0.200 |
| IBS | 107 | -0.253 | -0.259 |

Population-level tests:

- one-way ANOVA: `p = 1.75 × 10^-7`
- Kruskal-Wallis: `p = 3.33 × 10^-6`

These differences are not evidence that one of these populations has intrinsically higher or lower CAD incidence. PRS distributions can shift with allele frequencies, LD, discovery-GWAS composition, and calibration.

## chr22 pilot vs genome-wide PRS

For the same 503 individuals:

- Pearson `r = 0.0973`
- Spearman `rho = 0.1215`
- linear-model `R² = 0.00948`

The weak relationship is a useful internal lesson: a single-chromosome pilot can validate the workflow but captures very little of the variation in the complete autosomal score.
