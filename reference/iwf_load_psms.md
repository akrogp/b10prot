# Load PSMs from mzIdentML Files

This function loads Peptide-to-Spectrum Matches (PSMs) from mzIdentML
files in a specified directory, combines them, and adds a standard
column for the PSM score and protein reference.

## Usage

``` r
iwf_load_psms(path = ".", pattern = ".mzid", psm_score = NULL, verbose = FALSE)
```

## Arguments

- path:

  A character string specifying the directory where the mzIdentML files
  are located. Default is the current working directory (`"."`).

- pattern:

  A character string representing the file pattern to search for.
  Default is `".mzid"`, which targets mzIdentML files.

- psm_score:

  An optional character string specifying the column name to be used as
  the PSM score. If `NULL`, the last column in the loaded data will be
  used as the PSM score.

- verbose:

  Logical, whether to print progress messages (default TRUE).

## Value

A `data.frame` containing the combined PSM data, with added columns for
the PSM score (`psmScore`) and the protein reference (`proteinRef`).
