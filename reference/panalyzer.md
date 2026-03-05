# Run PAnalyzer on Peptide-to-Protein Data

Executes the PAnalyzer tool on a provided peptide-to-protein dataset.
This function writes the input data to a temporary file, runs the
PAnalyzer tool using Java, and reads the output into R.

## Usage

``` r
panalyzer(pep2prot)
```

## Arguments

- pep2prot:

  A data frame containing peptide-to-protein mappings. The data frame
  must include the columns `peptideRef` and `proteinRef`.

## Value

A data frame with the original data and additional columns:

- `peptideType`: Type of the peptide determined by PAnalyzer. Possible
  values include `"unique"`, `"discriminating"`, and
  `"non-discriminating"`.

- `proteinType`: Type of the protein determined by PAnalyzer. Possible
  values include `"conclusive"`, `"indistinguishable"`, `"ambiguous"`,
  and `"non-conclusive"`. All the proteins in a group will have the same
  `proteinType` value.

- `groupRef`: Group reference from PAnalyzer.

## See also

- [PAnalyzer: A software tool for protein inference in shotgun
  proteomics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-288)
  for more information on the tool.
