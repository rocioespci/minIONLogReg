\name{LASTparser_SplitForwardReverse}
\alias{LASTparser_SplitForwardReverse}
\title{LASTparser_SplitForwardReverse}
\usage{
  LASTparser_SplitForwardReverse(inFile, ref_label)
}
\description{
Function used for parsing .maf files. It splits .maf file into two files: one for forward reads (+) one for reverse reads (-).
}
\arguments{
  \item{inFile}{Character. Path to original file to be splitted. }
  \item{ref_label}{Character. Label used for refrence sequence in LAST alignment (it's the name indicated in every line where the reference sequence is.)}
}
\details{
This is a function used internally by parse_maf_to_alignment. It stores files in the directory indicated by inFile (working directory if none is specified).

Depends on awk.
}
\examples{


}
