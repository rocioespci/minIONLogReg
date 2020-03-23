\name{count_nucleotides_per_qscore}
\alias{count_nucleotides_per_qscore}
\title{count_nucleotides_per_qscore}
\usage{
count_nucleotides_per_qscore(inFile_prefix)
}
\description{
It counts in a mutiple sequence alignment the number of occurrences of each combination of position - nucleotide - Qscore.
}
\arguments{
\item{inFile_prefix}{Character. Name of the file to analyse, without the ".maf".}
}
\details{
It loads output from parse_maf_to_alignment.

It saves a file named inFile_prefix+"_countsPNQ.txt with a table with the following columns: position-nucleotide-quality-counts
}
\examples{
}

