\name{check_both_links}
\alias{check_both_links}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Verify if the provided QCTool's and Gene Cutter's results URL links are usable
}
\description{
This function verifies if the QCTool's and Gene Cutter's results URL links are usable.
}
\usage{
check_both_links(template_filename)
}
\arguments{
  \item{template_filename}{
the name of the Excel template file
}
}
\details{
This function is used at the beginning of the main function \preformatted{HIV_IntegrityAnalysis()}.
}
\value{
None
}
\references{
QCTool website: \url{https://www.hiv.lanl.gov/content/sequence/QC/index.html}
\cr Gene Cutter website: \url{https://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html}
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
