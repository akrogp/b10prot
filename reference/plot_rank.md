# Plot Rank of Decoy Scores

This function creates a rank plot of decoy scores based on various LP
(coLogarithm of Probability) metrics, including LPM, LPS, LPF, and LPG
scores (LPGM, LPGS, LPGF). The plot displays the ranked scores of decoys
with a reference line for comparison.

## Usage

``` r
plot_rank(data)
```

## Arguments

- data:

  A data frame containing identification data, including columns for
  decoy status (`isDecoy`) and any of the different LP metrics (LPM,
  LPS, LPF, LPGM, LPGS, LPGF).

## Value

A ggplot object showing the rank plot of decoy scores for the different
metrics. Each score type is displayed in a separate facet with the rank
plotted on the x-axis and the score on the y-axis. The red diagonal line
represents a reference for ideal ranking.
