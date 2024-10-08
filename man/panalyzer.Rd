% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/panalyzer.R
\name{panalyzer}
\alias{panalyzer}
\title{Run PAnalyzer on Peptide-to-Protein Data}
\usage{
panalyzer(pep2prot)
}
\arguments{
\item{pep2prot}{A data frame containing peptide-to-protein mappings. The
data frame must include the columns \code{peptideRef} and \code{proteinRef}.}
}
\value{
A data frame with the original data and additional columns:
\itemize{
\item \code{peptideType}: Type of the peptide determined by PAnalyzer. Possible values
include \code{"unique"}, \code{"discriminating"}, and \code{"non-discriminating"}.
\item \code{proteinType}: Type of the protein determined by PAnalyzer. Possible values
include \code{"conclusive"}, \code{"indistinguishable"}, \code{"ambiguous"}, and
\code{"non-conclusive"}. All the proteins in a group will have the same
\code{proteinType} value.
\item \code{groupRef}: Group reference from PAnalyzer.
}
}
\description{
Executes the PAnalyzer tool on a provided peptide-to-protein dataset. This
function writes the input data to a temporary file, runs the PAnalyzer tool
using Java, and reads the output into R.
}
\seealso{
\itemize{
\item \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-288}{PAnalyzer: A software tool for protein inference in shotgun proteomics} for more information on the tool.
}
}
