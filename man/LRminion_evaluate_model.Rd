\name{LRminion_evaluate_model}
\alias{LRminion_evaluate_model}
\title{LRminion_evaluate_model}
\usage{
LRminion_evaluate_model(filename_to_evaluate,filename_wildtype, LAST_reference_label, length_sequencing, ORF_ini,ORF_end,reference_sequence, auxiliary_files=NULL, verbose=T)
}
\description{
This is the main function to evaluate nanopore reads of a gene library.
}
\arguments{
\item{filename_to_evaluate}{Character vector. Files (.maf) to evaluate.}
\item{filename_wildtype}{Character vector of length one. Name of .maf with the reads of the reference sequence.}
  \item{LAST_reference_label}{Character. Label used for refrence sequence in LAST alignment (it's the name indicated in every line where the reference sequence is.)}
  \item{length_sequencing}{Numeric. Length of the region of interest in the sequences.}
  \item{ORF_ini}{Numeric. Position at the reference sequence where the region of interest in the sequence starts. }
  \item{ORF_end}{Numeric. Position at the reference sequence where the region of interest in the sequence ends.}
  \item{reference_sequence}{Character vector. Reference sequence in a vector (one nucleotide one element of the vector, in order.)}
  \item{auxiliary_files}{Should auxiliary files created in the process be removed or moved to another folder? NULL will remove them. A path will move them to that folder. }
  \item{verbose}{Logical. Should messages of the progress of the analysis be printed in the screen?}
}
\details{

}
\examples{

}

