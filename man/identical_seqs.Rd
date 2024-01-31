\name{identical_seqs}
\alias{identical_seqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count the number of different nucleotides between the aligned sequences
}
\description{
This function returns the number of different nucleotides between the aligned sequences.
}
\usage{
identical_seqs(matrix)
}
\arguments{
  \item{matrix}{
matrix of sequences aligned with MAFFT
}
}
\details{
Is used to identify the clones and potential clones of a sequence. This function is used in the main function \samp{Clonality_Analysis()}.
}
\value{
An integer representing the number of different nucleotides between to sequences
}
% \references{
% ProSeq-IT website: \url{https://psd.cancer.gov/tools/pvs_annot.php}
% }
\author{
Audrée Lemieux
}
\note{
%%  ~~further notes~~
}

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
