# Part 2 — From SNP selection to LD-aware models

This section moves from classical clumping + thresholding to methods that model LD and re-estimate SNP effects.

- **PRSice:** automates C+T threshold search.
- **LDpred2:** Bayesian LD-aware posterior effect estimation with chain diagnostics.
- **PRS-CS:** Bayesian regression with continuous shrinkage priors and an external LD reference.

The main conceptual transition is from asking **“which SNPs should I keep?”** to asking **“given LD, what posterior weight should each SNP receive?”**
