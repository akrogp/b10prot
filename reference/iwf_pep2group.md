# Assign Peptides to Protein Groups

This function processes peptide-to-protein group relationships by
summarizing the number of discriminating and total peptides per protein,
and then assigning peptides to their respective protein groups. Proteins
within a group are ordered by the number of discriminating peptides,
followed by the total number of peptides.

## Usage

``` r
iwf_pep2group(pep2prot2group, extra_pep_cols = c())
```

## Arguments

- pep2prot2group:

  A data frame containing peptide-to-protein-to-group relations. Must
  include at least `peptideRef`, `proteinRef`, `groupRef`,
  `peptideType`, and `proteinType`. It can also contain optional columns
  such as `shared`, `isDecoy`, `pepScore`, `LP`, or `qval`.

- extra_pep_cols:

  A character vector of additional peptide-level column names to include
  in the summarization step (if present in the input data). These are
  combined with a default set of optional columns (`shared`, `isDecoy`,
  `pepScore`, `LP`, `qval`). Columns not present in the data are
  silently ignored.

## Value

A data frame where peptides are assigned to protein groups. In addition
to the original peptide and group identifiers, the result contains:

- `proteinCount`: Number of proteins in each group.

- `proteinRefs`: Concatenation of all protein references in the group.

- `proteinMaster`: First protein of the group (after ordering).

- Any optional peptide-level columns present in the input (e.g.
  `shared`, `isDecoy`, `pepScore`, `LP`, `qval`).

## See also

[`iwf_grouping`](https://akrogp.github.io/b10prot/reference/iwf_grouping.md)
for generating peptide-to-protein-to-group relations, and
[`panalyzer`](https://akrogp.github.io/b10prot/reference/panalyzer.md)
for protein grouping.
