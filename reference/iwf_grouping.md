# Perform Protein Grouping Based on Peptide-to-Protein Relations

This function builds protein groups using the peptides that pass a
specified peptide-level FDR threshold. Peptides that do not pass the
threshold are grouped separately and assigned a negative group
identifier. It utilizes the PAnalyzer algorithm to infer protein groups
and assigns peptide and protein types.

## Usage

``` r
iwf_grouping(pep2prot, threshold = 0.01)
```

## Arguments

- pep2prot:

  A data frame containing peptide-to-protein relations, including the
  columns `peptideRef`, `proteinRef`, `qval`, and `isDecoy`.

- threshold:

  A numeric value specifying the peptide-level FDR threshold for
  grouping (default is 0.01).

## Value

A data frame with the inferred protein groups, with additional columns:

- `peptideType`: Type of peptide (`"unique"`, `"discriminating"`,
  `"non-discriminating"`).

- `proteinType`: Type of protein (`"conclusive"`, `"indistinguishable"`,
  `"ambiguous"`, `"non-conclusive"`).

- `groupRef`: Group identifier for the proteins.

- `shared`: The number of groups matched by each peptide.

## See also

[`panalyzer`](https://akrogp.github.io/b10prot/reference/panalyzer.md)
for protein grouping,
[`iwf_pep2level`](https://akrogp.github.io/b10prot/reference/iwf_pep2level.md)
for obtaining peptide-to-protein relations.
