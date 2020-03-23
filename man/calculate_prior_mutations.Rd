\name{calculate_prior_mutations}
\alias{calculate_prior_mutations}
\title{calculate_prior_mutations}
\usage{
calculate_prior_mutations(rates.matrix, mean.n.mut, wildtype_sequence, save=F, filename="tablePriorMutations.txt")
}
\description{

}
\arguments{
\item{rates.matrix}{Numeric matrix (dimension 4x5). Mutation rate in ePCR process, in order: A C G T -.}
\item{mean.n.mut}{Numeric (length 1). Mean number of mutations. }
\item{wildtype_sequence}{Character vector. Reference sequence in a vector (one nucleotide one element of the vector, in order.)}
\item{save}{Logical. Should results be stored?}
\item{filename}{Character (length 1). Name of the file where to store the results.}
}

\details{
Compute prior probabilities or having an actual mutation.
}


