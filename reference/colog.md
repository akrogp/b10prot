# Compute Logarithmic Transformation for Probabilities

This function computes a custom logarithmic transformation of
probabilities. If the probability is less than `LIMIT_PROB`, it returns
`LIMIT_COLOG`. Otherwise, it returns the negative base-10 logarithm of
the probability.

## Usage

``` r
colog(prob)
```

## Arguments

- prob:

  A numeric vector of probabilities.

## Value

A numeric vector where values less than `LIMIT_PROB` are replaced by
`LIMIT_COLOG`, and other values are transformed using `-log10()`.

## Examples

``` r
# Example with a single probability value
colog(0.001)
#> [1] 3

# Example with a vector of probabilities
colog(c(0.001, 0.05, 0.1))
#> [1] 3.00000 1.30103 1.00000
```
