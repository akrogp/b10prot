
<!-- README.md is generated from README.Rmd. Please edit that file -->

# b10prot <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->
<!-- badges: end -->

b10prot is an R package designed for the analysis of proteomics data,
specifically focusing on protein identification. It is developed as part
of the [EhuB10](https://ehubio.ehu.eus/) initiative, a collaborative
effort between several research groups from the University of the Basque
Country (UPV/EHU). The name b10prot is a reference to bioinformatics and
proteomics, with “b10” representing “bio” in a way that reflects both
biology and binary code.

This package is built with the aim of simplifying the integration of our
latest research into proteomics data analysis workflows. It works with
data in a “tidy” format, following principles similar to those of the
`tidyverse`.

## Key Features

- **Protein Inference** using the
  [PAnalyzer](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-288)
  algorithm:
  - `panalyzer` runs the PAnalyzer algorithm on peptide-to-protein data.
  - `plot_groups` plots PAnalyzer protein groups composition.
- **Rank Identifications** using the
  [LPGF](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00819) score:
  - `lpg` calculates the different LP Gamma (LPG) scores, including the
    recommended LPGF score.
  - `plot_rank` plots decoy scores vs their rank to check for an uniform
    distribution.
- **FDR Estimation** including the refined
  [FDRr](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00819)
  technique:
  - `target_decoy_approach` calculates p-values and q-values based on
    the traditional target-decoy approach.
  - `refined_fdr` computes different FDR estimations using a competitive
    approach between target and decoy identifications.

## Identification Workflow

The b10prot package includes a set of functions (with the `iwf_` prefix)
specifically designed to streamline the protein identification workflow.
These functions are designed to work with data in a “tidy” format,
following principles similar to those of the `tidyverse`. This means
that the data should be organized in a way that each type of observation
is stored in its own column and each row represents a single
observation.

This workflow is based on two main types of data:

- **Identification Lists** containing a list of identifications with
  their scores:
  - `iwf_load_psms` load PSMs from mzIdentML files.
  - `iwf_psm2pep` aggregates PSMs into peptides.
  - `lpg` collapses relationships into a list of identifications
    including LPG scores.
- **Identification Relationships** between lower-level (e.g., peptide)
  and higher-level (e.g., protein) identifications:
  - `iwf_pep2level` maps peptides to the specified level.
  - `iwf_grouping` performs protein grouping based on peptide-to-protein
    relations.
  - `iwf_pep2group` creates peptide-to-group relations from protein
    grouping relations.

## Installation

You can install the development version of b10prot from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("akrogp/b10prot")
```

## Tutorials

You can learn more in `vignette("b10prot")`.
