# Summarize PAnalyzer Results

Provides a summary of the PAnalyzer results by grouping the data based
on whether a protein is a target or decoy, and by its protein type. The
summary includes the number of distinct proteins and protein groups for
both target and decoy entries.

## Usage

``` r
# S3 method for class 'panalyzer'
summary(object, ...)
```

## Arguments

- object:

  A data frame containing the results from the `panalyzer` function.
  This data frame must contain the following columns:

  - `isDecoy`: Logical indicator of whether the protein is a decoy
    (`TRUE`) or a target (`FALSE`).

  - `proteinType`: The type of protein as classified by PAnalyzer (e.g.,
    `"conclusive"`, `"indistinguishable"`, `"ambiguous"`,
    `"non-conclusive"`).

  - `proteinRef`: Protein reference identifier.

  - `groupRef`: Protein group reference.

- ...:

  Additional arguments (currently ignored).

## Value

A data frame summarizing the PAnalyzer results with the following
columns:

- `Type`: The type of protein (`"conclusive"`, `"indistinguishable"`,
  `"ambiguous"`, `"non-conclusive"`).

- `TargetProteins`: The number of target proteins of the given type.

- `DecoyProteins`: The number of decoy proteins of the given type.

- `TargetGroups`: The number of target protein groups of the given type.

- `DecoyGroups`: The number of decoy protein groups of the given type.
