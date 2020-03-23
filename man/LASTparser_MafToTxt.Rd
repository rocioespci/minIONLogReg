
\name{LASTparser_MafToTxt}
\alias{LASTparser_MafToTxt}
\title{LASTparser_MafToTxt}
\usage{
LASTparser_MafToTxt(inFile_prefix, ref_label)
}
\description{
Function used for parsing .maf files. It splits .maf file into five:
reference sequence, reference sequence's initial position in alignment, aligned sequence, Qscore, name of read sequences.
}
\arguments{
\item{inFile_prefix}{Character. Name of the file to analyse, without the ".maf".}
\item{ref_label}{Character. Label used for refrence sequence in LAST alignment (it's the name indicated in every line where the reference sequence is).}
}
\details{
This is a function used internally by parse_maf_to_alignment. It stores files in the working directory.

Depends on awk.
}
\examples{


}

