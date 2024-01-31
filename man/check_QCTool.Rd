\name{check_QCTool}
\alias{check_QCTool}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Verify if the QCTool's summary text file provided is usable
}
\description{
This function verifies if the QCTool's summary text file provided is indeed a text file, if it can be found in your current directory, and if it contains the appropriate header.
}
\usage{
check_QCTool(QCTool_summary)
}
\arguments{
  \item{QCTool_summary}{
the name of the summary text (.txt) file downloaded from QCTool's output
}
}
\details{
This function is used at the beginning of the main function \preformatted{HIV_IntegrityAnalysis()}.
}
\value{
\preformatted{TRUE} if the QCTool's summary text file is usable, or \preformatted{FALSE} otherwise.
}
\references{
QCTool website: \url{https://www.hiv.lanl.gov/content/sequence/QC/index.html}
}
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
