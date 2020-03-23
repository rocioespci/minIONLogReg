\name{evaluate_fits}
\alias{evaluate_fits}
\title{evaluate_fits}
\usage{
  evaluate_fits(inFile_prefix, data_fits, reference_sequence, ascii)
}
\description{
Evaluates Qscores of new reads with the logistic regression built on the reference sequence.
}
\arguments{
  \item{inFile_prefix}{Character. Name of the file to analyse, without the ".maf".}
  \item{data_fits}{Data frame including all fits, as R object returned by fit_logistic_regression function.}
  \item{reference_sequence}{Character vector. Reference sequence in a vector (one nucleotide one element of the vector, in order.)}
  \item{ascii}{data frame with ascii translation from symbols to Qscores.}
}
\details{
Output: it saves a file named inFile_prefix+"_FrequenciesModelPriors.txt" with the new probabilities of being right.
}
\examples{

}

