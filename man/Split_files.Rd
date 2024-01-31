\name{Split_files}
\alias{Split_files}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Split the FASTA file
}
\description{
This function splits the FASTA file into shorter FASTA files, splitted by donor (one file per donor)
}
\usage{
Split_files(FASTA_file, donors)
}
\arguments{
  \item{FASTA_file}{
name of the FASTA file containing all the aligned sequences.
}
  \item{donors}{
vector string containing the names of the donors. Use \samp{c()} if more than one donor.
}
}
\details{
This function is used in the main function \samp{Clonality_Analysis()}.
}
\value{
None
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
