# Aggregate PSMs into Peptides

This function aggregates Peptide-to-Spectrum Matches (PSMs) into
peptides, selecting the best PSM score for each peptide and maintaining
decoy information. Peptides that are present in both the target and
decoy subsets are removed from the final results.

## Usage

``` r
iwf_psm2pep(psms, lower_better = TRUE)
```

## Arguments

- psms:

  A `data.frame` containing the PSM data, which must include columns
  `peptideRef`, `psmScore`, and `isDecoy`.

- lower_better:

  A logical value indicating whether a lower PSM score is better.
  Default is `TRUE`.

## Value

A `data.frame` containing one row per unique peptide with columns:

- `peptideRef`:

  The peptide identifier.

- `pepScore`:

  The best score for the peptide.

- `isDecoy`:

  A logical indicating whether the peptide belongs to the decoy set.

## See also

[iwf_load_psms](https://akrogp.github.io/b10prot/reference/iwf_load_psms.md)
for loading PSMs.
