% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoring.R
\name{global_fdr}
\alias{global_fdr}
\title{Compute Global False Discovery Rate (FDR)}
\usage{
global_fdr(data)
}
\arguments{
\item{data}{A data frame containing the identification data. It must include
a logical column \code{isDecoy} indicating whether each row corresponds to a decoy identification.}
}
\value{
A data frame with the total counts of target and decoy identifications,
and the global FDR expressed as a percentage.
}
\description{
This function calculates the global False Discovery Rate (FDR) as a percentage
based on the total number of target and decoy identifications in the data.
}
\examples{
# Example usage with a sample dataset
sample_data <- data.frame(
  isDecoy = c(TRUE, FALSE, FALSE, TRUE, FALSE)
)
global_fdr(sample_data)

}
