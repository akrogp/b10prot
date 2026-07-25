# Picked Protein Group FDR

Implements the picked protein group false discovery rate (FDR) strategy
as described by The et al. (2022).

## Usage

``` r
picked_gfdr(data, score, lower_better = TRUE, affix = "_REVERSED")
```

## Arguments

- data:

  A tibble containing protein group information. Must include at least
  the columns groupRef, proteinMaster, proteinRefs, isDecoy and the
  scoring column provided in score.

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

A tibble containing the subset of protein groups that passed the picked
protein group filtering and the subsequent target-decoy FDR analysis.
The output will typically include the original columns from data and any
additional FDR-related metrics added by target_decoy_approach().

## Details

Protein groups are first sorted by a scoring column, and for each
leading protein any subsequent group containing it within its members is
removed.

Internally, the function:

1.  Normalizes protein identifiers by removing the affix (decoy tag).

2.  Orders protein groups according to the chosen score and direction.

3.  Compares each group against subsequent groups to do the competition.

4.  Retains only the winner groups.

5.  Applies a target-decoy approach to estimate FDR.

## References

The M, Samaras P, Kuster B, Wilhelm M. Reanalysis of ProteomicsDB Using
an Accurate, Sensitive, and Scalable False Discovery Rate Estimation
Approach for Protein Groups. Mol Cell Proteomics. 2022
Dec;21(12):100437.
[doi:10.1016/j.mcpro.2022.100437](https://doi.org/10.1016/j.mcpro.2022.100437)

## See also

[`target_decoy_approach()`](https://akrogp.github.io/b10prot/reference/target_decoy_approach.md)

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
library(stringr)

df <- tibble(
  groupRef = c("a", "b", "c", "d", "e", "f_REVERSED", "f"),
  proteinMaster = c("a", "b", "c", "d", "e", "f_REVERSED", "f"),
  proteinRefs = c(
    "a;x",
    "a_REVERSED;b_REVERSED",
    "b_REVERSED;c_REVERSED;x_REVERSED",
    "d",
    "e;d",
    "f_REVERSED",
    "f"
  ),
  score = 7:1
) %>%
  mutate(isDecoy = str_detect(proteinRefs, "_REVERSED"))

result <- b10prot:::picked_gfdr(df, score, lower_better = FALSE)
glimpse(result)
#> Rows: 3
#> Columns: 12
#> $ groupRef      <chr> "a", "d", "f_REVERSED"
#> $ proteinMaster <chr> "a", "d", "f_REVERSED"
#> $ proteinRefs   <chr> "a;x", "d", "f_REVERSED"
#> $ score         <int> 7, 4, 2
#> $ isDecoy       <lgl> FALSE, FALSE, TRUE
#> $ decoys        <int> 0, 0, 1
#> $ targets       <int> 1, 2, 2
#> $ target        <int> 1, 2, 2
#> $ pval          <dbl> 0.5, 0.5, 0.5
#> $ LP            <dbl> 0.30103, 0.30103, 0.30103
#> $ FDR           <dbl> 0.0, 0.0, 0.5
#> $ qval          <dbl> 0.0, 0.0, 0.5
```
