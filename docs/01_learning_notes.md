# Learning notes — the project in plain language

## Part 1: learn what the score is doing

A PRS is mathematically simple: for each SNP, multiply the number of effect alleles carried by the SNP weight, then sum across SNPs.

What makes a real PRS difficult is everything around that equation.

The first exercises therefore isolate three failure modes:

- the effect allele can be wrong even when the software reports success;
- correlated SNPs can represent overlapping genetic information;
- the chosen P-value threshold changes which SNPs enter a classical C+T score.

## Part 2: stop treating SNPs as independent yes/no choices

PRSice automates the classical C+T strategy.

LDpred2 and PRS-CS introduce a different idea: use LD information to estimate posterior SNP effects. Instead of simply keeping or dropping a SNP, the model can shrink its effect continuously.

## Part 3: move to real CAD data

The real-data analysis was first restricted to chromosome 22 so that every step could be debugged cheaply.

After the pilot succeeded, the same logic was generalized to all autosomes. The final analysis produced genome-wide scores for 503 real 1000 Genomes EUR individuals.

The most important interpretation boundary is that these individuals are not a CAD case/control validation cohort. The analysis demonstrates genome-wide scoring and population-dependent score distributions, not clinical prediction accuracy.
