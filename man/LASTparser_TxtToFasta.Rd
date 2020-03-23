\name{LASTparser_TxtToFasta}
\alias{LASTparser_TxtToFasta}
\title{LASTparser_TxtToFasta}
\usage{
LASTparser_TxtToFasta(inFile_prefix, length_sequencing, ORF_ini, ORF_end)
}
\description{
Function used for parsing .maf files. It combines txt files produced by LASTparser_MafToTxt in an alignment stored at a .fasta file.

}
\arguments{
\item{inFile_prefix}{Character. Name of the file to analyse, without the ".maf".}
\item{length_sequencing}{Numeric. Length of the region of interest in the sequences.}
\item{ORF_ini}{Numeric. Position at the reference sequence where the region of interest in the sequence starts. }
\item{ORF_end}{Numeric. Position at the reference sequence where the region of interest in the sequence ends.}


}
\details{
This is a function used internally by parse_maf_to_alignment. It stores files in the working directory.

To each sequence:

- It constraints the analysis between ORF_ini and ORF_end.

- Fills with gaps ('-') the begining and end of the sequence (if there were no nucleotides read).

- Assigns deletions a Qscore equal to the average of the Qscores of nucleotides that are nearest neighbours.

- Writes two fasta files: One for the nucleotides, one for the Qscores.


}
\examples{

}

