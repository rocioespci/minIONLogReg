\name{parse_maf_to_alignment}
\alias{parse_maf_to_alignment}
\title{parse_maf_to_alignment}
\usage{
parse_maf_to_alignment(filename, LAST_reference_label, length_sequencing, ORF_ini,ORF_end, auxiliary_files=NULL)
}
\description{
This function parses .maf files into two .fasta files (one for nucleotides, one for Qscores).

Take into account that:

- I will store data in the working path or path indicated in the "inFile" input.

- It uses awk command.

}
\arguments{
  \item{inFile}{Character. Path to original file to be splitted. }
  \item{ref_label}{Character. Label used for refrence sequence in LAST alignment (it's the name indicated in every line where the reference sequence is.)}
  \item{length_sequencing}{Numeric. Length of the region of interest in the sequences.}
  \item{ORF_ini}{Numeric. Position at the reference sequence where the region of interest in the sequence starts. }
  \item{ORF_end}{Numeric. Position at the reference sequence where the region of interest in the sequence ends.}
  \item{auxiliary_files}{Should auxiliary files created in the process be removed or moved to another folder? NULL will remove them. A path will move them to that folder. }

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



