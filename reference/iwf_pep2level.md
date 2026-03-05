# Map Peptides to a Specified Reference Level

This function maps peptides to a specified reference level (e.g.,
protein, gene, or protein group) based on Peptide-to-Spectrum Matches
(PSMs). It also calculates how many reference entities (e.g., proteins,
genes, or groups) are matched by each peptide.

## Usage

``` r
iwf_pep2level(psms, levelRef)
```

## Arguments

- psms:

  A `data.frame` containing PSM data, including `peptideRef` (peptide
  identifier) and the reference level specified by `levelRef`.

- levelRef:

  The column name in `psms` that contains the reference level to which
  peptides should be mapped.

## Value

A `data.frame` containing:

- `peptideRef`:

  The peptide identifier.

- `levelRef`:

  The reference level identifier.

- `shared`:

  The number of reference entities matched by each peptide.
