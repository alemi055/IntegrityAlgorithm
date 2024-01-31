\name{check_integer}
\alias{check_integer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Verify if the \samp{analyzes} and \samp{threshold} values provided are integers
}
\description{
This function verifies if the \samp{analyzes} and \samp{threshold} values provided are integers
}
\usage{
check_integer(analyzes, threshold)
}
\arguments{
  \item{analyzes}{
integer. The functions to run. \samp{1}: QCTool only; \samp{2}: GeneCutter and ProSeq-IT; \samp{3}: IntegrateInfo only; \samp{4}: All
}
  \item{threshold}{
integer. The threshold of different nucleotides between two sequences to consider them as \strong{potential clones}.
}
}
\details{
This function is used at the beginning of the main functions \samp{HIV_IntegrityAnalysis()} and \samp{Clonality_Analysis()}.
}
\value{
None
}
% \references{
% QCTool website: \url{https://www.hiv.lanl.gov/content/sequence/QC/index.html}
% Gene Cutter website: \url{https://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html}
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
