# Refined Group-Level False Discovery Rate (FDRr) Calculation

This function computes refined False Discovery Rate (FDR) estimates at
the protein group level using a competitive approach. It extends
`refined_fdr` by handling multiple proteins per identification in the
target-decoy competition.

## Usage

``` r
refined_gfdr(data, score, lower_better = TRUE, affix = "_REVERSED")
```

## Arguments

- data:

  A data frame containing identification data. It must include:

  - `proteinMaster`: the main protein identifier for the group.

  - `proteinRefs`: all proteins associated with the identification,
    separated by ";".

  - `isDecoy`: logical indicating whether the identification is a decoy.

  - A score column used for ranking identifications.

- score:

  The column name of the score used to rank the identifications. This
  should be an unquoted column name.

- lower_better:

  Logical; if TRUE, lower scores indicate better identifications
  (default `TRUE`).

- affix:

  String indicating the suffix/prefix used to identify decoy entries
  (default `"_REVERSED"`).

## Value

A data frame containing the original data along with additional columns:

- FDRn:

  Normal FDR estimation as cumulative minimum (q-value).

- FDRp:

  Picked FDR estimation as cumulative minimum (q-value).

- FDRr:

  Refined FDR estimation as cumulative minimum (q-value).

- to:

  Target-only identifications count at the protein group level.

- do:

  Decoy-only identifications count at the protein group level.

- td:

  Count of identifications with equal target and decoy scores.

- tb:

  Target-best identifications count at the protein group level.

- db:

  Decoy-best identifications count at the protein group level.

## Details

The function first expands multiple protein references per
identification into individual rows, removes decoy affixes, and then
determines which protein references are in competition with the
proteinMaster. It calculates the FDR estimates by summing target and
decoy identifications across regions (to, do, td, tb, db), following the
competitive approach described in the cited paper.

## See also

- [Refined FDR for Single
  Proteins](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00819)

- [Competitive Protein Group
  FDR](https://pubmed.ncbi.nlm.nih.gov/36328188)
