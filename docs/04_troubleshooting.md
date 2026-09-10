# Troubleshooting

Several issues came up while moving from toy examples to a real genome-wide PRS workflow. The points below are kept because they are directly useful for understanding PRS implementation and reproducibility.

## 1. Allele direction

An early toy exercise deliberately used an inconsistent effect allele without flipping the beta sign.

The scoring command still completed successfully, which made the lesson very clear:

> a technically valid score file does not guarantee biologically correct allele orientation.

This is why allele harmonization must be checked explicitly before scoring.

## 2. PRS-CS and recent NumPy behavior

The downloaded PRS-CS source encountered scalar-conversion errors under a recent NumPy environment.

The issue arose because one-element NumPy arrays were being used in places that expected scalar values. Minimal local compatibility edits converted those one-element arrays to scalar values before scalar-only calculations and output formatting.

The modified third-party PRS-CS source is not redistributed in this repository.

## 3. PRS-CS-auto `phi`

Passing the literal string:

```text
--phi=None
```

did not invoke PRS-CS-auto. PRS-CS attempted to parse `None` as a numeric value and failed.

PRS-CS-auto was correctly invoked by **omitting the `--phi` argument entirely**.

## 4. 1000 Genomes v5b variant IDs

During the chromosome-22 pilot, some 1000 Genomes Phase 3 v5b VCF records had a missing ID field (`.`).

Variant IDs therefore had to be recovered from chromosome, position, and allele information before the final scoring step.

This was one reason the later genome-wide implementation used a cleaner input route.

## 5. EBI bulk-download speed

The first genome-wide download attempt through EBI was very slow and uneven.

A public AWS 1000 Genomes Phase 3 v5a route was tested and was substantially faster from the available network environment. The successful genome-wide workflow therefore used the AWS v5a files consistently across chromosomes 1–22.

The superseded partial v5b download attempt is not included in this repository.
