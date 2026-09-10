# 08 — Real-world CAD chromosome-22 pilot

The pilot uses real CAD GWAS summary statistics, the PRS-CS 1000G EUR HapMap3 LD reference, and 503 European-ancestry 1000 Genomes individuals.

The purpose of the chromosome-22 analysis was **workflow validation**, not construction of a chromosome-22 CAD predictor.

Validated steps:

GWAS QC → HapMap3 position/allele harmonization → frequency QC → target genotype preparation → PRS-CS-auto → posterior weights → PLINK2 scoring → standardized individual PRS.

The final chr22 scoring set contained **16,431 SNPs** in **503 individuals**.
