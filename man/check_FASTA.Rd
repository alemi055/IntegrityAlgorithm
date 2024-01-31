\name{check_FASTA}
\alias{check_FASTA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Verify if the FASTA file provided is usable
}
\description{
This function verifies if the FASTA file provided, containing the aligned sequences, is indeed a FASTA file and if it can be found in your current directory.
}
\usage{
check_FASTA(FASTA_file)
}
\arguments{
  \item{ProseqIT_rx}{
name of the FASTA file containing all the aligned sequences
}
}
\details{
This function is used at the beginning of the main function \preformatted{HIV_IntegrityAnalysis()}.
}
\value{
\preformatted{TRUE} if the FASTA file is usable, or \preformatted{FALSE} otherwise.
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
