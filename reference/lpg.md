# LP Gamma (LPG) Metrics Calculation

This function calculates various LP Gamma (LPG) metrics for a specified
reference level, such as a protein, a gene or a protein group. The
metrics are derived from the coLogarithm of Probability (LP) of their
constituting peptides and include probabilities derived from maximum LP,
sum of LP values, and filtered LP sums using a specified peptide-level
FDR threshold.

## Usage

``` r
lpg(data, levelRef, threshold = 0.01, extra_cols = c())
```

## Arguments

- data:

  A data frame containing identification data, including peptide-related
  columns for coLogarithm of Probability (`LP`), q-values (`qval`), and
  a logical column `isDecoy`.

- levelRef:

  The column name of the reference level to group by, such as a protein
  or a gene identifier. This should be an unquoted column name.

- threshold:

  A numeric value representing the FDR threshold for peptide-level
  q-values (default is `0.01`).

- extra_cols:

  A character vector of optional columns to include in the summarisation
  step (if present in the input data). By default, it attempts to
  include `proteinType`, `proteinCount`, `proteinRefs`, and
  `proteinMaster`.

## Value

A data frame containing the calculated LPG metrics with the following
columns:

- isDecoy:

  Indicates whether the group contains any decoy identification.

- n:

  The total number of peptide identifications for the group.

- m:

  The number of peptide identifications with a q-value below the
  threshold.

- LPM:

  The maximum coLogarithm of Probability (`LP`) for the group.

- LPS:

  The sum of coLogarithm of Probability (`LP`) for the group.

- LPF:

  The sum of coLogarithm of Probability for identifications with a
  q-value below the threshold.

- LPGM:

  The LP Gamma value based on the maximum `LP`.

- LPGS:

  The LP Gamma value based on the sum of `LP` values.

- LPGF:

  The LP Gamma value based on the filtered sum of `LP` values for
  confident identifications.

## See also

- [Protein Probability Model for High-Throughput Protein Identification
  by Mass Spectrometry-Based
  Proteomics](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00819)
  for more information on the LPG scores.

## Examples

``` r
# Example usage with a sample dataset
sample_data <- data.frame(
  levelRef = c("P1", "P1", "P2", "P2", "P3"),
  LP = c(1.5, 2.0, 0.5, 1.0, 1.2),
  qval = c(0.01, 0.02, 0.005, 0.03, 0.01),
  isDecoy = c(FALSE, FALSE, TRUE, FALSE, FALSE)
)
lpg(sample_data, levelRef, threshold = 0.01)
#> # A tibble: 3 × 10
#>   levelRef isDecoy     n     m   LPM   LPS   LPF  LPGM  LPGS  LPGF
#>   <chr>    <lgl>   <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 P1       FALSE       2     1   2     3.5   1.5 1.70  2.54  1.20 
#> 2 P2       TRUE        2     1   1     1.5   0.5 0.721 0.851 0.199
#> 3 P3       FALSE       1     1   1.2   1.2   1.2 1.2   1.2   1.2  
```
