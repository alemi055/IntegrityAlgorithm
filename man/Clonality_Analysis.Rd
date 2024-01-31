\name{Clonality_Analysis}
\alias{Clonality_Analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assess the clonality of HIV proviruses
}
\description{
This function assesses the clonality of HIV proviruses. For donor and for each unique sequence, it returns the number and names of its clones and the number and names of its "potential clones", if applicable. It outputs the results in the "Clonality" folder, in the "FINAL_OUTPUT" folder).
}
\usage{
Clonality_Analysis(FASTA_file, donors, threshold = 5)
}
\arguments{
  \item{FASTA_file}{
name of the FASTA file containing all the aligned sequences.
}
  \item{donors}{
vector string containing the names of the donors. Use c() if more than one donor.
}
  \item{threshold}{
integer. The threshold of different nucleotides between two sequences to consider them as "potential clones". Default is 5
}
}
\details{
To analyze the data, the function first splits the FASTA file into shorter FASTA files, splitted by donor (one file per donor). It then creates a matrix of sequence length and paired-comparisons are only done with sequences with the same length, or \eqn{\302\261} the threshold. It realigns these sequences with MAFFT. If the number of different nucleotides between two sequences is 0, and their sequence length are identical, the second sequence is classified as a "clone". Instead, it is \eqn{\302\261} the threshold, it is classified as a "potential clone". That can later be confirmed using Geneious (see \url{https://github.com/alemi055/IntegrityAlgorithm#to-confirm-clones-and-potential-clones-in-geneious}).
}
\value{
None
}
% \references{
% QCTool website: \url{https://www.hiv.lanl.gov/content/sequence/QC/index.html}
% Gene Cutter website: \url{https://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html}
% ProSeq-IT website: \url{https://psd.cancer.gov/tools/pvs_annot.php}
% }
\author{
Audrée Lemieux
}
% \note{
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
%
% ## The function is currently defined as
% function (x)
% {
%   }
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory (show via RShowDoc("KEYWORDS")):
% % \keyword{ ~kwd1 }
% % \keyword{ ~kwd2 }
% % Use only one keyword per line.
% % For non-standard keywords, use \concept instead of \keyword:
% % \concept{ ~cpt1 }
% % \concept{ ~cpt2 }
% % Use only one concept per line.
