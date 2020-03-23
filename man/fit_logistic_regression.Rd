\name{fit_logistic_regression}
\alias{fit_logistic_regression}
\title{fit_logistic_regression}
\usage{
fit_logistic_regression(inFile_prefix, reference_sequence, prior_error, prior_mutation)
}
\description{
Fits logistic regression position and nucleotide dependent, for Qscore correction.
}
\arguments{
\item{inFile_prefix}{Character. Name of the file to analyse, without the ".maf".}
\item{reference_sequence}{Character vector. Reference sequence in a vector (one nucleotide one element of the vector, in order.)}
\item{prior_error}{Data frame with the columns position,nucleotide,prior.error as returned by calculate_prior_errors function. }
\item{prior_mutation}{Data frame with the columns  wt.base,nucleotide,p_mutation as returned by calculate_prior_mutations function.}

}
\details{

Take into account that it takes a long time to run.

It stores two Rdata files, one with the fittings (named inFile_prefix+"_apriori_fits.Rdata") and one with the fitted data (named inFile_prefix+"_apriori_data.Rdata").
}
\examples{
}

