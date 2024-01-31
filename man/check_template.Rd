\name{check_template}
\alias{check_template}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Verify if the Excel template file provided is usable
}
\description{
This function verifies if the Excel template file provided is indeed an Excel file, if it can be found in your current directory, and if it contains the appropriate tabs.
}
\usage{
check_template(template_filename)
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
\preformatted{TRUE} if the Excel template file is usable, or \preformatted{FALSE} otherwise.
}
% \references{
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
