# Compute Global False Discovery Rate (FDR)

This function calculates the global False Discovery Rate (FDR) as a
percentage based on the total number of target and decoy identifications
in the data.

## Usage

``` r
global_fdr(data)
```

## Arguments

- data:

  A data frame containing the identification data. It must include a
  logical column `isDecoy` indicating whether each row corresponds to a
  decoy identification.

## Value

A data frame with the total counts of target and decoy identifications,
and the global FDR expressed as a percentage.

## Examples

``` r
# Example usage with a sample dataset
sample_data <- data.frame(
  isDecoy = c(TRUE, FALSE, FALSE, TRUE, FALSE)
)
global_fdr(sample_data)
#> # A tibble: 1 × 3
#>   Target Decoy `Global FDR (%)`
#>    <int> <int>            <dbl>
#> 1      3     2             66.7
```
