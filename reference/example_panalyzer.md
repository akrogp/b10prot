# PAnalyzer Example Data

A dataset containing the PAnalyzer example presented in the original
paper.

## Usage

``` r
example_panalyzer
```

## Format

A data frame with 18 rows and the following 5 columns:

- peptideRef:

  Character vector of peptide references.

- proteinRef:

  Character vector of protein references.

- peptideType:

  Character vector indicating the type of peptide.

- proteinType:

  Character vector indicating the type of protein.

- groupRef:

  Integer vector of protein group references.

## Source

Data extracted from [PAnalyzer: A software tool for protein inference in
shotgun
proteomics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-288)
for demonstration purposes.

## Examples

``` r
data(example_panalyzer, package = "b10prot")
plot_groups(example_panalyzer, groupRefs = 1:5)

```
