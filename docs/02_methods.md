# Methods

## Data

### CAD GWAS

The real-world case study used CAD GWAS summary statistics corresponding to **GWAS Catalog accession GCST90132314** (Aragam et al., 2022) in a GRCh37-compatible representation.

### LD reference

PRS-CS used the **1000 Genomes European HapMap3 LD reference**.

### Target genotypes

The target genotype panel contained **503 European-ancestry individuals** from 1000 Genomes Phase 3, spanning CEU, FIN, GBR, IBS, and TSI.

## GWAS QC and harmonization

The final genome-wide pipeline applied the same logic chromosome by chromosome:

- autosomes 1–22;
- A/C/G/T biallelic SNPs;
- finite beta and standard error;
- `SE > 0`;
- GWAS MAF calculated as `min(EAF, 1-EAF)`;
- `MAF >= 0.01`;
- valid positive case and control counts;
- per-variant effective sample size:
  `N_eff = 4 / (1/N_case + 1/N_control)`;
- match to HapMap3 by chromosome and GRCh37 position;
- direct or swapped allele matches retained;
- beta sign flipped for swapped matches;
- absolute GWAS/reference MAF difference `<= 0.10`;
- duplicate matched rsIDs reduced to one retained record.

Across autosomes, the successful run retained **1,119,312 SNPs**.

## Target genotype preparation

Public 1000 Genomes Phase 3 autosomal genotype data were used as the target panel source.

The VCFs were subset to the predefined 503-person EUR list and harmonized positions, converted with PLINK2, assigned/recovered variant IDs as required, and reduced to the final chromosome-specific SNP sets.

Target allele QC required every retained SNP to be a direct or swapped allele match.

## PRS-CS-auto

The final workflow ran PRS-CS separately for chromosomes 1–22.

Because the CAD GWAS contains variant-specific sample sizes while PRS-CS expects a single `n_gwas`, the project used the **median post-QC effective sample size across the genome: 601,956** as a pragmatic approximation.

PRS-CS-auto was run chromosome by chromosome using the same harmonized SNP sets and ancestry-matched LD reference.

These MCMC settings were used to validate the workflow; a formal research analysis should evaluate convergence and sensitivity more thoroughly.

## Individual scoring

Posterior SNP effects were scored with PLINK2. Chromosome-specific `SCORE1_SUM` values were summed across chr1–22 to create the genome-wide PRS, then standardized to a Z-score across the 503 target individuals.

## Population analysis

Descriptive summaries were generated for CEU, FIN, GBR, IBS, and TSI. One-way ANOVA and Kruskal-Wallis tests were used to test for distributional differences among the five groups.

These tests evaluate **score distributions**, not CAD prevalence or clinical risk.
