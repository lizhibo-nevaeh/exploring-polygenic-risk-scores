# Limitations

1. **No phenotype validation.**  
   The 503 1000 Genomes EUR individuals are a genotype target panel, not an independent CAD case/control cohort. This project therefore does not estimate AUC, odds ratios per PRS SD, calibration, or clinical utility.

2. **Single `n_gwas` approximation.**  
   The CAD meta-analysis contains variant-specific sample sizes, while PRS-CS accepts one `n_gwas`. The genome-wide median effective sample size after QC was used as a pragmatic approximation.

3. **Demonstration-length MCMC.**  
   The final PRS-CS-auto run used 1000 iterations with 500 burn-in and thinning of 5. More formal analyses should examine convergence and sensitivity under longer chains.

4. **Ancestry dependence.**  
   PRS values depend on ancestry, LD reference choice, allele frequencies, discovery-GWAS composition, and calibration. Differences among CEU/FIN/GBR/IBS/TSI should not be read as direct differences in CAD incidence.

5. **Portfolio scope.**  
   This repository emphasizes learning, QC logic, and reproducible workflow design. It is not a clinical software package or a novel PRS method.
