\name{QCTool_analyzes}
\alias{QCTool_analyzes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyze the results from QCTool
}
\description{
This function analyzes the results from QCTool.

It first identifies the stop codons, then hypermutations (presence of hypermutation: \samp{Possible} or absence: \samp{Cannot determine}). This function creates a temporary file named \samp{Analyzed_QCTool.csv} in the \strong{tmp} folder.
}
\usage{
QCTool_analyzes(filename, summary, RefSeq = TRUE)
}
\arguments{
  \item{filename}{
the name of the Excel template file
}
  \item{summary}{
the name of the summary text (.txt) file downloaded from QCTool's output
}
  \item{RefSeq}{
logical. If \samp{TRUE}, the reference sequence is included in QCTool's and Gene Cutter's results. Default is \samp{TRUE}
}
}
\details{
To analyze the data, the function assesses the summary text (.txt) file downloaded from QCTool's output, as well as the URL link of the results (to access the stop codons).

This function is used in the main function \samp{HIV_IntegrityAnalysis()}.
}
\value{
None
}
\references{
QCTool website: \url{https://www.hiv.lanl.gov/content/sequence/QC/index.html}
}
\author{
Audrée Lemieux
}
\note{
\strong{Troubleshooting}

~~~~~~~~~~~~~~~~~~~~

\code{Warning in file(file, ifelse(append, "a", "w")): cannot open file 'tmp/Analyzed_QCTool.csv': No such file or directory}\cr
\code{Error in file(file, ifelse(append, "a", "w")): cannot open the connection}

This error occurs when the user does not have administrative rights. To fix this error, manually create the folders: \strong{tmp} and \strong{FINAL_OUTPUT}.
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
