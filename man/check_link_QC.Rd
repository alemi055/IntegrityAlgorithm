\name{check_link_QC}
\alias{check_link_QC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Verify if the provided QCTool's results URL link is usable
}
\description{
This function verifies if the QCTool's results URL link is a \samp{https://} link, is indeed from LANL's HIV Sequence Database, is indeed results from QCTool (and not another LANL's tool, for instance), and if the link has not expired.
}
\usage{
check_link_QC(hyperlink)
}
\arguments{
  \item{hyperlink}{
URL link of LANL's HIV Sequence Database's QCTool's results
}
}
\details{
This function is used at in the function \samp{check_both_links()}.
}
\value{
\samp{TRUE} if the QCTool's results URL link provided is usable, or \samp{FALSE} otherwise.
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
