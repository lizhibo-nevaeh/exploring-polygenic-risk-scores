# 02 — Genome-wide CAD PRS

This stage extends the validated chromosome-22 pilot to autosomes 1–22.

The scientific workflow is:

1. split the CAD GWAS by chromosome;
2. apply SNP-level QC and HapMap3 harmonization;
3. summarize the post-QC effective sample size;
4. prepare ancestry-matched target genotypes;
5. estimate posterior SNP effects with PRS-CS-auto;
6. calculate chromosome-specific scores with PLINK2;
7. sum chromosome scores into a genome-wide PRS;
8. standardize the scores and compare their distributions.

The completed analysis retained **1,119,312 HapMap3-matched SNPs** across autosomes and used a median effective GWAS sample size of **601,956** as the single `n_gwas` approximation required by PRS-CS.

Only the scientific analysis scripts and compact results needed to understand the workflow are included here. Compute-environment and job-scheduling files are intentionally excluded from the portfolio.
