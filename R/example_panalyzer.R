#' PAnalyzer Example Data
#'
#' A dataset containing the PAnalyzer example presented in the original paper.
#'
#' @format A data frame with 18 rows and the following 5 columns:
#' \describe{
#'   \item{peptideRef}{Character vector of peptide references.}
#'   \item{proteinRef}{Character vector of protein references.}
#'   \item{peptideType}{Character vector indicating the type of peptide.}
#'   \item{proteinType}{Character vector indicating the type of protein.}
#'   \item{groupRef}{Integer vector of protein group references.}
#' }
#'
#' @examples
#' data(example_panalyzer, package = "b10prot")
#' plot_groups(example_panalyzer, groupRefs = 1:5)
#'
#' @source Data extracted from [PAnalyzer: A software tool for protein inference in shotgun proteomics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-288) for demonstration purposes.
"example_panalyzer"
