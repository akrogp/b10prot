# Check for required columns in a data frame

This helper function verifies that all required columns exist in a data
frame or tibble. If any are missing, it throws an error.

## Usage

``` r
check_required_cols(data, required)
```

## Arguments

- data:

  A data frame or tibble.

- required:

  A character vector of required column names.

## Value

The same tibble if all required columns exist.
