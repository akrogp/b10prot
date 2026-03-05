# Plot Protein Groups Composition

Visualizes the composition of selected protein groups based on their
associated peptides. The function filters the specified protein groups,
arranges them, and displays a plot showing the peptides associated with
each protein in the group. Peptides are colored according to their
classification as unique, discriminating, non-discriminating, or
non-significant.

## Usage

``` r
plot_groups(panalyzer, groupRefs)
```

## Arguments

- panalyzer:

  A data frame containing the results from the `panalyzer` function. The
  data frame must include columns such as `groupRef`, `peptideRef`,
  `peptideType`, `proteinRef` and `proteinType`.

- groupRefs:

  A vector of group reference identifiers (`groupRef`) representing the
  protein groups to be visualized.

## Value

A ggplot object visualizing the selected protein groups and their
associated peptides.

## See also

- [PAnalyzer: A software tool for protein inference in shotgun
  proteomics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-288)
  for an example plot.

## Examples

``` r
data(example_panalyzer, package = "b10prot")
plot_groups(example_panalyzer, groupRefs = 1:5)

```
