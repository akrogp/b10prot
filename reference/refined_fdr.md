# Refined False Discovery Rate (FDR) Calculation

This function computes refined False Discovery Rate (FDR) estimates
using a competitive approach between target and decoy identifications.
It provides three types of FDR calculations: FDRn, FDRp, and FDRr, which
adjust for different competitive scenarios between targets and decoys.

## Usage

``` r
refined_fdr(data, levelRef, score, lower_better = TRUE, affix = "_REVERSED")
```

## Arguments

- data:

  A data frame containing the identification data, including columns for
  the reference level, a score, and whether each identification is a
  decoy (`isDecoy`).

- levelRef:

  The column name containing the reference level for each identification
  (e.g., protein or gene reference). This should be an unquoted column
  name.

- score:

  The column name of the score used to rank the identifications. This
  should be an unquoted column name.

- lower_better:

  A logical value indicating whether lower scores are better (default is
  `TRUE`).

- affix:

  A string indicating the suffix of prefix used to identify decoy
  entries in the reference level column. Default is `"_REVERSED"`.

## Value

A data frame with the original data and additional columns for the
refined FDR estimates:

- FDRn:

  Normal FDR estimation as cumulative minimum (q-value).

- FDRp:

  Picked FDR estimation as cumulative minimum (q-value).

- FDRr:

  Refined FDR estimationas cumulative minimum (q-value).

- to:

  Target-only identifications count.

- do:

  Decoy-only identifications count.

- td:

  Count of identifications with the same target and decoy scores.

- tb:

  Target-best identifications count.

- db:

  Decoy-best identifications count.

## See also

- [Protein Probability Model for High-Throughput Protein Identification
  by Mass Spectrometry-Based
  Proteomics](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00819)
  for more information on the refined FDR estimation.

## Examples

``` r
# Example usage with a sample dataset
sample_data <- data.frame(
  proteinRef = c("P1", "P1_REVERSED", "P2", "P3", "P3_REVERSED"),
  score = c(0.1, 0.2, 0.3, 0.5, 0.4),
  isDecoy = c(FALSE, TRUE, FALSE, FALSE, TRUE)
)
refined_fdr(sample_data, levelRef = proteinRef, score = score, lower_better = TRUE)
#> # A tibble: 5 × 11
#>   proteinRef  score isDecoy    to    do    td    tb    db  FDRn  FDRp  FDRr
#>   <chr>       <dbl> <lgl>   <int> <int> <int> <int> <int> <dbl> <dbl> <dbl>
#> 1 P1            0.1 FALSE       1     0     0     0     0 0       0   0    
#> 2 P1_REVERSED   0.2 TRUE        0     0     0     1     0 0.5     0   0    
#> 3 P2            0.3 FALSE       1     0     0     1     0 0.5     0   0    
#> 4 P3_REVERSED   0.4 TRUE        1     1     0     1     0 0.667   0.5 0.5  
#> 5 P3            0.5 FALSE       1     0     0     1     1 0.667   0.5 0.667
```
