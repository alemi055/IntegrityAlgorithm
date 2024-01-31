\name{check_link_GC}
\alias{check_link_GC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Verify if the provided Gene Cutter's results URL link is usable
}
\description{
This function verifies if the Gene Cutter's results URL link is a \'https://\' link, is indeed from LANL's HIV Sequence Database, is indeed results from Gene Cutter (and not another LANL's tool, for instance), and if the link has not expired.
}
\usage{
check_link_GC(hyperlink)
}
\arguments{
  \item{hyperlink}{
URL link of LANL's HIV Sequence Database's Gene Cutter's results
}
}
\details{
This function is used at in the function check_both_links().
}
\value{
TRUE if the Gene Cutter's results URL link provided is usable, or FALSE otherwise.
}
\references{
Gene Cutter website: \url{https://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html}
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
